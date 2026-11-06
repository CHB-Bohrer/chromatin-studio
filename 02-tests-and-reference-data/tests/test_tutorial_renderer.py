"""Narration controls the visual clock and cannot be cut to fit old scenes."""

import copy
from contextlib import redirect_stderr, redirect_stdout
import importlib.util
import io
import json
import math
from pathlib import Path
import struct
import sys
import tempfile
import unittest
from unittest.mock import patch
import wave


RENDERER_PATH = Path(__file__).resolve().parents[1] / "tutorial" / "render_tutorial.py"
SPEC = importlib.util.spec_from_file_location("tutorial_renderer_under_test", RENDERER_PATH)
renderer = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(renderer)


class TutorialRendererTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.captures = self.root / "captures"
        self.captures.mkdir()
        self.data = {"normal_duration_seconds": 8.5, "scenes": [
            self.scene("long-speech", 2.5), self.scene("short-speech", 6),
        ]}
        self.manifest = {"engine": "test PCM fixture", "scenes": []}
        self.pcm = {}
        for scene, frame_count in zip(self.data["scenes"], (168_731, 53_461)):
            path = self.root / (scene["id"] + ".wav")
            # Distinct final samples make truncation detectable, even if the
            # output has subsequently been padded to the expected duration.
            payload = b"".join(struct.pack("<hh", (i % 801) - 400, 400 - (i % 801))
                               for i in range(frame_count))
            self.pcm[scene["id"]] = payload
            with wave.open(str(path), "wb") as audio:
                audio.setnchannels(2)
                audio.setsampwidth(2)
                audio.setframerate(48_000)
                audio.writeframes(payload)
            self.manifest["scenes"].append({
                "id": scene["id"], "file": path.name, "source_text": scene["narration"],
                "wav_sha256": renderer.sha256(path), "duration_seconds": frame_count / 48_000,
                "voice": "test voice", "words_per_minute": 130,
                "time_stretch": False, "speech_truncated": False,
            })
            for keyframe in scene["keyframes"]:
                (self.captures / keyframe["capture"]).write_bytes(b"fixture capture")
        self.manifest_path = self.root / "narration-manifest.json"
        self.save_manifest()

    @staticmethod
    def scene(identifier, duration):
        return {"id": identifier, "title": identifier.title(),
                "narration": "Complete spoken text for " + identifier,
                "normal_duration_seconds": duration,
                "keyframes": [{"capture": identifier + "-a.png", "duration": duration * .4},
                              {"capture": identifier + "-b.png", "duration": duration * .6}],
                "captions": [{"start": 0, "end": duration * .3, "text": "First key point"},
                             {"start": duration * .35, "end": duration, "text": "Final key point"}]}

    def save_manifest(self):
        self.manifest_path.write_text(json.dumps(self.manifest))

    def timeline(self, fps=30, pause=1.25):
        return renderer.narration_timeline(self.data, self.manifest_path, fps, pause)

    def test_long_narration_extends_scene_without_mutating_source_storyboard(self):
        original = copy.deepcopy(self.data)
        adjusted = self.timeline()
        scene = adjusted["scenes"][0]
        measured = self.manifest["scenes"][0]["duration_seconds"]
        expected = math.ceil((measured + 1.25) * 30) / 30
        self.assertGreater(measured, original["scenes"][0]["normal_duration_seconds"])
        self.assertEqual(scene["normal_duration_seconds"], expected)
        self.assertGreaterEqual(scene["narration_audio"]["reading_pause_seconds"], 1.25)
        self.assertLess(scene["narration_audio"]["reading_pause_seconds"], 1.25 + 1 / 30)
        self.assertEqual(scene["source_duration_seconds"], 2.5)
        self.assertEqual(scene["narration_audio"]["acceleration"], 1)
        self.assertFalse(scene["narration_audio"]["speech_truncated"])
        self.assertFalse(scene["narration_audio"]["time_stretch"])
        self.assertEqual(adjusted["playback_speed"], 1)
        self.assertEqual(self.data, original)

    def test_keyframes_captions_and_exported_chapters_share_the_complete_timeline(self):
        adjusted = self.timeline(fps=24, pause=.7)
        total = sum(scene["normal_duration_seconds"] for scene in adjusted["scenes"])
        self.assertEqual(adjusted["normal_duration_seconds"], total)
        self.assertEqual(adjusted["rendered_duration_seconds"], total)
        expected_caption_times = []
        offset = 0
        for source, scene in zip(self.data["scenes"], adjusted["scenes"]):
            duration = scene["normal_duration_seconds"]
            self.assertAlmostEqual(duration * 24, round(duration * 24))
            self.assertAlmostEqual(sum(k["duration"] for k in scene["keyframes"]), duration)
            self.assertAlmostEqual(scene["keyframes"][0]["duration"] / duration, .4)
            intervals = renderer.scene_intervals(scene)
            self.assertAlmostEqual(intervals[0][0], 0)
            self.assertAlmostEqual(intervals[-1][1], duration)
            self.assertAlmostEqual(sum(end - start for start, end, _, _ in intervals), duration)
            for previous, following in zip(intervals, intervals[1:]):
                self.assertAlmostEqual(previous[1], following[0])
            factor = duration / source["normal_duration_seconds"]
            for before, caption in zip(source["captions"], scene["captions"]):
                self.assertGreaterEqual(caption["start"], 0)
                self.assertLessEqual(caption["end"], duration + 1e-9)
                self.assertLess(caption["start"], caption["end"])
                self.assertAlmostEqual(caption["start"], before["start"] * factor)
                self.assertAlmostEqual(caption["end"], before["end"] * factor)
                expected_caption_times.append((offset + caption["start"], offset + caption["end"]))
            offset += duration
        chapters, metadata = renderer.write_accessibility(adjusted, self.root, "tutorial", speed=1)
        self.assertEqual(chapters[0]["start"], 0)
        self.assertAlmostEqual(chapters[-1]["end"], total)
        self.assertEqual(chapters[0]["end"], chapters[1]["start"])
        self.assertEqual(json.loads((self.root / "tutorial-chapters.json").read_text()), chapters)
        self.assertIn(f"END={round(total * 1000)}", metadata.read_text())
        vtt = (self.root / "tutorial.vtt").read_text()
        srt = (self.root / "tutorial.srt").read_text()
        for start, end in expected_caption_times:
            self.assertIn(renderer.timestamp(start) + " --> " + renderer.timestamp(end), vtt)
            self.assertIn(renderer.timestamp(start, True) + " --> " + renderer.timestamp(end, True), srt)

    def test_stale_missing_or_modified_narration_is_rejected(self):
        self.data["scenes"][0]["narration"] += " An added sentence."
        with self.assertRaisesRegex(ValueError, "stale"):
            self.timeline()
        self.data["scenes"][0]["narration"] = self.manifest["scenes"][0]["source_text"]
        removed = self.manifest["scenes"].pop()
        self.save_manifest()
        with self.assertRaisesRegex(ValueError, "missing"):
            self.timeline()
        self.manifest["scenes"].append(removed)
        self.save_manifest()
        path = self.root / self.manifest["scenes"][0]["file"]
        path.write_bytes(path.read_bytes()[:-4])
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            self.timeline()

    def test_corrupt_wav_header_is_rejected_even_with_matching_hash(self):
        record = self.manifest["scenes"][0]
        path = self.root / record["file"]
        path.write_bytes(b"not a valid WAV header")
        record["wav_sha256"] = renderer.sha256(path)
        self.save_manifest()
        with self.assertRaises((ValueError, wave.Error, EOFError)):
            self.timeline()

    def test_truncated_pcm_payload_is_rejected_even_with_matching_hash(self):
        record = self.manifest["scenes"][0]
        path = self.root / record["file"]
        # Keep the WAV header's full duration, but remove half the speech data.
        path.write_bytes(path.read_bytes()[:-len(self.pcm[record["id"]]) // 2])
        record["wav_sha256"] = renderer.sha256(path)
        self.save_manifest()
        with self.assertRaises(ValueError):
            self.timeline()

    def test_manifest_cannot_claim_stretched_or_truncated_speech_is_complete(self):
        for flag in ("time_stretch", "speech_truncated"):
            for value in (True, None):
                with self.subTest(flag=flag, value=value):
                    self.manifest["scenes"][0][flag] = value
                    self.save_manifest()
                    with self.assertRaisesRegex(ValueError, "unaccelerated and complete"):
                        self.timeline()
                    self.manifest["scenes"][0][flag] = False

    def test_prepare_audio_preserves_every_speech_sample_and_only_appends_silence(self):
        scene = self.timeline()["scenes"][0]
        work = self.root / "audio-work"
        work.mkdir()
        output, speech = renderer.prepare_audio(scene, work, renderer.imageio_ffmpeg.get_ffmpeg_exe())
        with wave.open(str(output), "rb") as audio:
            self.assertEqual((audio.getnchannels(), audio.getsampwidth(), audio.getframerate()), (2, 2, 48_000))
            self.assertEqual(audio.getnframes(), round(scene["normal_duration_seconds"] * 48_000))
            delivered = audio.readframes(audio.getnframes())
        source = self.pcm[scene["id"]]
        self.assertEqual(delivered[:len(source)], source)
        self.assertEqual(set(delivered[len(source):]), {0})
        self.assertEqual(speech, scene["narration_audio"])

    def test_prepare_audio_refuses_a_visual_duration_shorter_than_speech(self):
        scene = self.timeline()["scenes"][0]
        scene["normal_duration_seconds"] = scene["narration_audio"]["duration_seconds"] - .01
        with patch.object(renderer, "run") as subprocess_run:
            with self.assertRaisesRegex(ValueError, "refusing to truncate"):
                renderer.prepare_audio(scene, self.root, "ffmpeg")
            subprocess_run.assert_not_called()

    def test_narrated_cli_refuses_acceleration_before_rendering(self):
        for speed in ("2", "1.25", "0.75"):
            with self.subTest(speed=speed):
                error = io.StringIO()
                with patch.object(sys, "argv", [str(RENDERER_PATH), "--speed", speed, "--validate-only"]), \
                        patch.object(renderer, "run") as subprocess_run, redirect_stderr(error):
                    with self.assertRaises(SystemExit) as raised:
                        renderer.main()
                self.assertEqual(raised.exception.code, 2)
                self.assertIn("recorded pace", error.getvalue())
                subprocess_run.assert_not_called()

    def test_silent_validation_retains_explicit_caption_only_speed_option(self):
        storyboard = self.root / "storyboard.json"
        storyboard.write_text(json.dumps(self.data))
        output = io.StringIO()
        with patch.object(sys, "argv", [str(RENDERER_PATH), "--silent", "--speed", "2",
                                        "--validate-only", "--storyboard", str(storyboard),
                                        "--captures", str(self.captures)]), redirect_stdout(output):
            self.assertEqual(renderer.main(), 0)
        report = json.loads(output.getvalue())
        self.assertTrue(report["valid"])
        self.assertEqual(report["delivered_seconds"], self.data["normal_duration_seconds"] / 2)


if __name__ == "__main__":
    unittest.main()

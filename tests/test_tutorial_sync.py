"""Each tutorial chapter must stay self-consistent so one chapter can be edited alone.

These checks read files only. They do not need the voice model or the renderer.
A failure names the chapter to fix; see tutorial/README.md, "Update one chapter".
"""
import json
import unittest
from pathlib import Path

TUTORIAL = Path(__file__).resolve().parents[1] / "tutorial"


class TutorialSyncTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.storyboard = json.loads((TUTORIAL / "storyboard.json").read_text(encoding="utf-8"))
        cls.scenes = cls.storyboard["scenes"]

    def test_scene_ids_are_unique(self):
        # Ids are stable names, not positions: the list order is the playback order.
        identifiers = [scene["id"] for scene in self.scenes]
        self.assertEqual(len(set(identifiers)), len(identifiers))

    def test_transcript_is_exactly_the_storyboard_narration(self):
        expected = "\n\n".join(scene["title"] + "\n" + scene["narration"] for scene in self.scenes)
        actual = (TUTORIAL / "transcript.txt").read_text(encoding="utf-8")
        blocks = [block.strip() for block in actual.strip().split("\n\n")]
        self.assertEqual(len(blocks), len(self.scenes), "transcript.txt needs one title+narration block per scene.")
        for scene, block in zip(self.scenes, blocks):
            self.assertEqual(block, scene["title"] + "\n" + scene["narration"],
                             f"transcript.txt differs from storyboard.json for scene {scene['id']}.")
        self.assertEqual(actual.strip(), expected)

    def test_every_referenced_capture_exists(self):
        for scene in self.scenes:
            names = {scene["primary_capture"]} | {frame["capture"] for frame in scene["keyframes"]}
            for name in names:
                self.assertTrue((TUTORIAL / "captures" / name).is_file(), f"{scene['id']}: missing capture {name}")

    def test_captions_and_keyframes_fit_the_reference_duration(self):
        for scene in self.scenes:
            duration = scene["normal_duration_seconds"]
            captions = scene["captions"]
            self.assertEqual(captions[0]["start"], 0, scene["id"])
            for earlier, later in zip(captions, captions[1:]):
                self.assertEqual(earlier["end"], later["start"], f"{scene['id']}: caption gap or overlap")
            self.assertEqual(captions[-1]["end"], duration, scene["id"])
            self.assertAlmostEqual(sum(frame["duration"] for frame in scene["keyframes"]), duration,
                                   msg=f"{scene['id']}: keyframe durations must sum to the reference duration")

    def test_recorded_narration_matches_the_current_text(self):
        manifest = json.loads((TUTORIAL / "narration-neural" / "narration-manifest.json").read_text(encoding="utf-8"))
        recorded = [entry["source_text"] for entry in manifest["scenes"]]
        stale = [scene["id"] for scene, text in zip(self.scenes, recorded) if scene["narration"] != text]
        self.assertEqual(len(recorded), len(self.scenes), "Narration audio is missing for at least one chapter.")
        self.assertEqual(stale, [], "Narration audio is stale for these chapters; rerun narrate_neural.py --scene ID.")


if __name__ == "__main__":
    unittest.main()

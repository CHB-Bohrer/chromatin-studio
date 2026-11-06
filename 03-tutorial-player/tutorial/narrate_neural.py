#!/usr/bin/env python3
"""Generate local Kokoro narration at a measured natural pace, without atempo.

Install requirements-narration.txt in a separate environment. Use --download
once to fetch pinned, checksum-verified public assets (~327 MB); subsequent
runs are offline. Output is one WAV and metadata file per scene, plus a manifest.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import re
import time
import urllib.request

import numpy as np
import onnxruntime as ort
import soundfile as sf
from kokoro_onnx import Kokoro

HERE = Path(__file__).resolve().parent
REPOSITORY = "onnx-community/Kokoro-82M-v1.0-ONNX"
REVISION = "1939ad2a8e416c0acfeecc08a694d14ef25f2231"
ASSETS = {
    "onnx/model.onnx": "8fbea51ea711f2af382e88c833d9e288c6dc82ce5e98421ea61c058ce21a34cb",
    "voices/af_heart.bin": "d583ccff3cdca2f7fae535cb998ac07e9fcb90f09737b9a41fa2734ec44a8f0b",
    "voices/af_bella.bin": "f69d836209b78eb8c66e75e3cda491e26ea838a3674257e9d4e5703cbaf55c8b",
}

def sha256(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()

def prepare_assets(root, download=False):
    root = Path(root)
    for name, expected in ASSETS.items():
        path = root / name
        if not path.is_file():
            if not download:
                raise ValueError(f"Missing {path}. Run once with --download to fetch public model assets.")
            path.parent.mkdir(parents=True, exist_ok=True)
            temporary = path.with_suffix(path.suffix + ".part")
            url = f"https://huggingface.co/{REPOSITORY}/resolve/{REVISION}/{name}"
            print("Downloading " + name, flush=True)
            with urllib.request.urlopen(url, timeout=120) as response, temporary.open("wb") as stream:
                while chunk := response.read(1024 * 1024):
                    stream.write(chunk)
            if sha256(temporary) != expected:
                raise ValueError("Downloaded asset failed its SHA-256 check: " + name)
            temporary.replace(path)
        if sha256(path) != expected:
            raise ValueError("Model asset failed its SHA-256 check: " + name)
    # Raw voice tensors are numeric float32 data. No pickle deserialization.
    voices = {name: np.fromfile(root / "voices" / (name + ".bin"), dtype=np.float32).reshape(-1, 1, 256)
              for name in ("af_heart", "af_bella")}
    packed = root / "voices-selected.npz"
    np.savez(packed, **voices)
    return root / "onnx/model.onnx", packed

def spoken_text(text):
    """Expand scientific initialisms for speech; preserve the source transcript."""
    substitutions = {
        "TFs": "transcription factors", "TF": "transcription factor",
        "IDR": "I D R", "DNA": "D N A", "RNA": "R N A",
        "Kd": "K D", "GPU": "G P U", "CPU": "C P U", "MSD": "M S D",
    }
    for word, replacement in substitutions.items():
        text = re.sub(r"\b" + re.escape(word) + r"\b", replacement, text)
    # Explicit separation prevents the English front end dropping the vowel.
    text = re.sub(r"\bcoactivator(s?)\b", r"co-activator\1", text, flags=re.IGNORECASE)
    return text

def count_words(text):
    return len(re.findall(r"\b[\w]+(?:[-'][\w]+)*\b", text))

def narrate(engine, source, voice, target_wpm, initial_speed=.90):
    text = spoken_text(source)
    count = count_words(text)
    if not count:
        raise ValueError("Narration is empty.")
    speed = initial_speed
    attempts = []
    candidates = []
    phonemes = engine.tokenizer.phonemize(text, "en-us")
    overrides = []
    if re.search(r"\bcohesin\b", text, flags=re.IGNORECASE):
        original, corrected = "kˈoʊhɪsˌɪn", "koʊhˈiːsɪn"
        occurrences = phonemes.count(original)
        if occurrences:
            phonemes = phonemes.replace(original, corrected)
            overrides.append({"word": "cohesin", "phonemes": corrected, "occurrences": occurrences})
    # Change the model's duration control, never stretch or speed up the waveform.
    for attempt in range(5):
        started = time.monotonic()
        audio, sample_rate = engine.create(phonemes, voice=voice, speed=speed, lang="en-us", is_phonemes=True,
                                          sentence_pause=.4, clause_pause=.12)
        audio = np.asarray(audio, dtype=np.float32)
        if len(audio) < sample_rate / 2 or not np.isfinite(audio).all():
            raise ValueError("Neural speech returned empty or non-finite audio.")
        duration = len(audio) / sample_rate
        wpm = count * 60 / duration
        attempts.append({"model_speed": speed, "seconds": duration, "words_per_minute": wpm,
                         "inference_seconds": time.monotonic() - started})
        candidates.append((abs(wpm - target_wpm), audio, sample_rate, duration, wpm, speed))
        if abs(wpm - target_wpm) <= 7 or attempt == 4:
            break
        next_speed = min(1.35, max(.65, speed * target_wpm / wpm))
        slower = [a["model_speed"] for a in attempts if a["words_per_minute"] < target_wpm]
        faster = [a["model_speed"] for a in attempts if a["words_per_minute"] > target_wpm]
        if slower and faster and max(slower) < min(faster):
            next_speed = (max(slower) + min(faster)) / 2
        if abs(next_speed - speed) < .005:
            break
        speed = next_speed
    _, audio, sample_rate, duration, wpm, speed = min(candidates, key=lambda candidate: candidate[0])
    peak = float(np.max(np.abs(audio)))
    if peak <= .001:
        raise ValueError("Neural speech is unexpectedly silent.")
    # Only attenuate clipping risk; do not compress dynamic range or alter timing.
    gain = min(1., .98 / peak)
    audio = audio * gain
    pad = np.zeros(round(.25 * sample_rate), dtype=np.float32)
    audio = np.concatenate([pad, audio, pad])
    info = {"voice": voice, "source_text": source, "spoken_text": text, "word_count": count,
            "target_words_per_minute": target_wpm, "words_per_minute": wpm,
            "model_speed": speed, "duration_seconds": len(audio) / sample_rate,
            "speech_seconds": duration, "sample_rate": sample_rate, "padding_seconds": .5,
            "peak_before_attenuation": peak, "gain": gain, "attempts": attempts,
            "within_target_pace": abs(wpm - target_wpm) <= 10,
            "phoneme_overrides": overrides,
            "time_stretch": False, "speech_truncated": False}
    return audio, sample_rate, info

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-dir", type=Path, default=Path("work/voice-models"))
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--storyboard", type=Path, default=HERE / "storyboard.json")
    parser.add_argument("--output-dir", type=Path, default=HERE / "narration-neural")
    parser.add_argument("--scene", action="append")
    parser.add_argument("--voice", choices=("af_heart", "af_bella"), default="af_heart")
    parser.add_argument("--target-wpm", type=float, default=155.)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if not math.isfinite(args.target_wpm) or not 120 <= args.target_wpm <= 185:
        parser.error("Use a target pace between 120 and 185 words per minute.")
    if not 1 <= args.threads <= 16:
        parser.error("Use 1–16 inference threads.")
    data = json.loads(args.storyboard.read_text())
    scenes = data["scenes"]
    ids = {s["id"] for s in scenes}
    if args.scene and set(args.scene) - ids:
        parser.error("Unknown scene ID.")
    selected = [s for s in scenes if not args.scene or s["id"] in args.scene]
    if any(not re.fullmatch(r"[a-z0-9][a-z0-9-]*", s["id"]) for s in selected):
        parser.error("Scene IDs must be filename-safe.")
    model, voices = prepare_assets(args.model_dir, args.download)
    options = ort.SessionOptions()
    options.intra_op_num_threads = args.threads
    options.inter_op_num_threads = 1
    session = ort.InferenceSession(str(model), sess_options=options, providers=["CPUExecutionProvider"])
    engine = Kokoro.from_session(session, str(voices))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    previous_path = args.output_dir / "narration-manifest.json"
    previous = json.loads(previous_path.read_text()) if previous_path.is_file() else {}
    records = []
    versions = {name: importlib.metadata.version(name) for name in
                ("kokoro-onnx", "onnxruntime", "numpy", "soundfile", "phonemizer", "espeakng-loader")}
    for scene in selected:
        path = args.output_dir / (scene["id"] + ".wav")
        metadata = path.with_suffix(".json")
        fingerprint = hashlib.sha256(json.dumps({"source": scene["narration"], "voice": args.voice,
            "target_wpm": args.target_wpm, "script": sha256(__file__), "versions": versions,
            "model_revision": REVISION}, sort_keys=True).encode()).hexdigest()
        if not args.force and path.is_file() and metadata.is_file():
            cached = json.loads(metadata.read_text())
            if cached.get("fingerprint") == fingerprint and cached.get("wav_sha256") == sha256(path):
                records.append(cached)
                print("Reusing " + scene["id"], flush=True)
                continue
        print("Narrating " + scene["id"] + " with " + args.voice, flush=True)
        audio, rate, info = narrate(engine, scene["narration"], args.voice, args.target_wpm)
        sf.write(str(path), audio, rate, subtype="PCM_16")
        info.update({"id": scene["id"], "file": path.name, "fingerprint": fingerprint,
                     "wav_sha256": sha256(path), "generator_script_sha256": sha256(__file__)})
        metadata.write_text(json.dumps(info, indent=2) + "\n")
        records.append(info)
        print(f"  {info['duration_seconds']:.2f} seconds; {info['words_per_minute']:.1f} words/min", flush=True)
    # Updating a scene retains previously verified audio for other unchanged scenes.
    # Each row keeps its own generator hash, so mixed revisions remain explicit.
    preserved = []
    if args.scene and previous.get("model_revision") == REVISION and previous.get("versions") == versions:
        by_id = {s["id"]: s for s in scenes}
        for old in previous.get("scenes", []):
            if old["id"] in args.scene or old["id"] not in by_id:
                continue
            path = args.output_dir / (old["id"] + ".wav")
            if (old.get("source_text") == by_id[old["id"]]["narration"] and old.get("voice") == args.voice
                    and old.get("target_words_per_minute") == args.target_wpm and path.is_file()
                    and old.get("wav_sha256") == sha256(path)):
                old.setdefault("generator_script_sha256", previous.get("narration_script_sha256"))
                records.append(old)
                preserved.append(old["id"])
        order = {s["id"]: index for index, s in enumerate(scenes)}
        records.sort(key=lambda r: order[r["id"]])
    manifest = {"schema_version": 1, "engine": "Kokoro-82M v1.0 / ONNX Runtime CPU",
                "voice": args.voice, "model_repository": REPOSITORY, "model_revision": REVISION,
                "model_assets_sha256": ASSETS, "versions": versions, "time_stretch": False,
                "narration_script_sha256": sha256(__file__), "storyboard_sha256": sha256(args.storyboard),
                "preserved_scene_ids": preserved, "complete": len(records) == len(scenes),
                "total_seconds": sum(r["duration_seconds"] for r in records), "scenes": records,
                "sources": ["https://huggingface.co/hexgrad/Kokoro-82M",
                            "https://huggingface.co/onnx-community/Kokoro-82M-v1.0-ONNX",
                            "https://github.com/thewh1teagle/kokoro-onnx"]}
    (args.output_dir / "narration-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Generated {len(records)} scenes, {manifest['total_seconds']:.2f} seconds total.", flush=True)

if __name__ == "__main__":
    main()

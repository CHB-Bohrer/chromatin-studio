# Editable Chromatin Studio tutorial

The tutorial uses actual application captures, on-screen key points and natural-paced neural narration. It is a guided UI walkthrough, not a continuous trajectory recording. The default player runs at 1× and offers a slower pace or chapter navigation.

The voice is synthesized locally with Kokoro, using the `af_heart` voice and a target of 155 words per minute. No audio acceleration is applied. Each chapter lasts as long as its measured narration, followed by a short reading pause. The previous fast `-2x` export remains available as a legacy file; the player now uses `chromatin-studio-tutorial.mp4`.

## Editable sources

- **storyboard.json** — chapter order, titles, narration, key points, capture actions and relative keyframe timings.
- **captures/** — actual UI screenshots with stable filenames.
- **transcript.txt** — the complete narration text; keep it synchronized with script edits.
- **narrate_neural.py** — local voice generation, pronunciation expansion and audio caching.
- **render_tutorial.py** — scene timing, video rendering, player and media packaging.
- **narration-neural/** — reusable scene WAVs and their source-text, voice, pace and checksum records.

Storyboard times are reference timings. The renderer scales each scene's visual intervals to the measured narration plus a 1.25-second pause, rounding the chapter duration up to a video frame. It never squeezes speech into a fixed scene length. The final timeline is saved in **render-timeline.json**; chapter and timed-text files use those final times. The on-screen notes and optional text track summarize each step; the full spoken script is available separately in the transcript.

## Generate narration

Keep narration dependencies in a separate Python 3.12 environment. Install **requirements-narration.txt** there. On this workspace, that environment is `../../work/voice-env` relative to the project directory; existing model assets are in `../../work/voice-models`.

From the project directory:

```sh
../../work/voice-env/bin/python tutorial/narrate_neural.py \
  --model-dir ../../work/voice-models \
  --output-dir tutorial/narration-neural
```

For a fresh setup, add `--download` once to fetch approximately 327 MB of public model and voice assets. They are pinned to a repository revision and checked against SHA-256 hashes. Subsequent runs work offline. No account or API key is required. The model sources and versions are recorded in the narration manifest.

Optional settings include `--voice af_bella`, `--target-wpm 145`, and `--threads 4`. The model's duration control sets speaking pace; the waveform is not sped up afterward. Unchanged audio is reused automatically. Use `--force` only when a full rebuild is intended.

## Render

Use a Python environment with **Pillow** and **imageio-ffmpeg**, such as this workspace's `../../work/sim-env/bin/python`.
On Windows, set `PYTHONUTF8=1` for both narration and rendering; the renderer falls back to `C:WindowsFontsrial.ttf`. After rendering, make sure `render-manifest.json` and `render-timeline.json` hold repo-relative `narration-neural/...` paths, not absolute local ones.

```sh
python tutorial/render_tutorial.py \
  --work-dir /path/to/work/tutorial-render \
  --capture-provenance tutorial/capture-provenance.json
```

The default narration manifest is `tutorial/narration-neural/narration-manifest.json`. Supply `--narration-manifest /path/to/narration-manifest.json` to use another directory. The renderer verifies scene text, audio hashes, complete PCM data and measured durations before accepting it. Stale or missing narration must be regenerated.

Other options:

```sh
python tutorial/render_tutorial.py --validate-only
python tutorial/render_tutorial.py --width 1280 --height 720
python tutorial/render_tutorial.py --reading-pause 2
python tutorial/render_tutorial.py --silent
```

The silent version uses reference timings without narration. Narrated exports require speed 1; viewers can change playback speed in the player. This prevents accidental double acceleration.

Generated files include the MP4, timed key points in SRT/VTT, chapter JSON/FFmetadata, poster, render manifest, resolved timeline and offline **index.html** player with **player.js**. Keep the player files beside the video and captions. No external browser libraries are required.

## Update one chapter

1. Edit its wording or capture actions in **storyboard.json**.
2. Recapture only affected UI views, retaining their filenames.
3. Update the matching narration in **transcript.txt**.
4. Run the narration command with `--scene 07-regulation` (using the edited scene ID). Other verified scenes stay in the complete manifest.
5. Render a preview with `python tutorial/render_tutorial.py --scene 07-regulation`.
6. Run the full render command to assemble the player, updated chapter times and video.

Use the same work directory and render settings to reuse cached clips. Renderer changes invalidate scene caches. Changed chapter lengths can also invalidate later clips because their progress indicator changes.

A keyframe can specify `focus_rect: [x, y, width, height]` in source-image pixels. Keep important status labels visible. Without a crop, the full screenshot is letterboxed. PNG and JPEG content are supported, including captures with a `.png` filename containing JPEG pixels.

## Keeping the tutorial current (checklist for future AI contributors)

Run this whenever a feature, preset, limit or default changes.

1. **Single source of truth.** `storyboard.json` holds narration, captions and capture actions; `transcript.txt` must contain the same narration, in storyboard scene order (nothing checks this automatically, so compare them by hand). Edit both together.
2. **Scale vocabulary.** A bead is a coarse genomic interval whose size is a setting (`chromatin.bp_per_bead`). At 200 bp it is a coarse nucleosome unit (occupied/vacant). At 10 kb it is a bin of many nucleosomes. Never say "each bead is X" without naming the preset; the opening overview uses the 20 Mb preset (10 kb), while chapters 02–17 use 200 bp presets. Check against `chromatin_studio/schema.py` `PRESETS` and `examples/*.json`.
3. **Claims audit.** Compare every number and limit (bead counts, bounds, GPU scope, calibration status) with `AGENTS.md`, `docs/PROJECT_STATUS.md` and `CONTRACT.md`. Narration must not claim calibration, validation or full-GPU execution that the status file does not.
4. **Mark audio stale.** Changing a scene's text invalidates its `narration-neural/<scene>.wav`; regenerate with `narrate_neural.py --scene <id>`, then re-render. If the voice environment is unavailable, say so in `PROJECT_STATUS.md` and the `status.json` tutorial note rather than leaving the video described as current.
5. **Recapture** only the UI views whose labels or layout changed, keeping filenames; record provenance.
6. **Update handoff.** Add a dated entry to `docs/PROJECT_STATUS.md` and keep `status.json` consistent.

## Capture and scientific review

Use an isolated recording instance or preserve and restore the current run. Only actual engine steps should create demonstrated bindings, contacts and RNA. Distinguish three recording models: the coarse 20 Mb overview, the 80 kb manuscript-hypothesis nucleosome-unit example, and the focused 25.6 kb nucleosome-access example. Chapters **03-genome** and **16-remodelers** explicitly switch to the 25.6 kb example; **04-proteins** explicitly returns to the 80 kb example. Preserve separate model/run provenance for these captures. The storyboard now contains **17 chapters** (chapter 17 covers the whole nucleus and tethering). The v0.10 set of 16 added a live 1D binding walkthrough and a dedicated accessibility/remodeler chapter to confinement, the live regulatory rate map, body responses, within-interval promoter timing and isolated manuscript studies. The rendered duration comes from measured narration; do not reuse the previous video duration after rebuilding.

Verify the distinctions between draft and applied models, learning illustrations and measurements, affinity and observed occupancy, basal transcription and TF-dependent gains, coarse IDRs and molecular detail, synthetic calibration and empirical evidence, implementation tests and published-figure replication, and cluster export and completed GPU execution.

For the new 1D chapter, capture **Live binding**, then use **Binding detail** at the 200 bp mapping to show individual sites and bound molecules. Select a site for the **Association / t** and **Release / t** table. The remodeler chapter needs actual views of the initial vacant interval, species capture/release factors, eviction/reassembly/sliding controls, and **Recorded binding & remodeling events**. Only engine-generated events are acceptable; do not fabricate a slide or release if none occurred. The recent-event buffer is bounded, and turnover timestamps resolve reaction intervals rather than every microscopic transition. Preserve these labels in the capture.

Stable new/updated capture names are **03-genome.png**, **03-genome-detail.png**, **03-genome-affinities.png**, **16-nucleosome-gaps.png**, **16-access-rates.png**, **16-remodelers.png** and **16-remodeler-events.png**. Recapture **04-proteins.png** to reflect the current capture/release editor labels; unchanged chapters can reuse verified captures. Keep raw captures unchanged after recording.

Review the final video at 1× after changing the script. Technical pronunciation is editable in `spoken_text()` within **narrate_neural.py**. Inspect audio, key-point readability, transitions and chapter seeking. Source hashes document which files were rendered; they do not independently establish scientific validity.

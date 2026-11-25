#!/usr/bin/env python3
"""Render actual-UI captures with natural-paced, pre-generated neural narration.

Narrated scenes follow measured WAV durations, with a short reading pause.
Audio is never accelerated or truncated. Requires Pillow and imageio-ffmpeg.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import html
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import wave
from PIL import Image, ImageDraw, ImageFont, ImageOps
import imageio_ffmpeg

HERE = Path(__file__).resolve().parent
FONT_ROOT = Path("/System/Library/Fonts/Supplemental")
C = {"bg":"#080b0f","panel":"#121923","line":"#2b3847","text":"#eef3f7",
     "muted":"#aab8c7","blue":"#79b7de","coral":"#ee7687"}

def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

def run(command):
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode:
        raise RuntimeError("Command failed: "+Path(command[0]).name+"\n"+result.stderr[-6000:])

def load_storyboard(path):
    data = json.loads(Path(path).read_text())
    scenes = data.get("scenes")
    if not isinstance(scenes,list) or not scenes:
        raise ValueError("Storyboard needs a nonempty scenes list.")
    ids = set()
    for scene in scenes:
        name = scene["id"]
        if not re.fullmatch(r"[a-z0-9][a-z0-9-]*",name) or name in ids:
            raise ValueError("Scene IDs must be unique filename-safe identifiers.")
        ids.add(name)
        duration = float(scene["normal_duration_seconds"])
        if not math.isfinite(duration) or duration <= 0:
            raise ValueError(name+": duration must be positive.")
        if not math.isclose(sum(float(k["duration"]) for k in scene["keyframes"]),duration,abs_tol=1e-7):
            raise ValueError(name+": keyframes do not sum to scene duration.")
        previous = 0
        for caption in scene["captions"]:
            start,end = float(caption["start"]),float(caption["end"])
            if start < previous or end <= start or end > duration:
                raise ValueError(name+": captions must be ordered and inside the scene.")
            previous = end
        for keyframe in scene["keyframes"]:
            filename = Path(keyframe["capture"])
            if filename.is_absolute() or ".." in filename.parts or float(keyframe["duration"])<=0:
                raise ValueError(name+": invalid capture path or duration.")
    total = sum(float(s["normal_duration_seconds"]) for s in scenes)
    if not math.isclose(total,float(data.get("normal_duration_seconds",total)),abs_tol=1e-7):
        raise ValueError("Declared duration does not match the scenes.")
    return data

def font(size,bold=False):
    candidates=[FONT_ROOT/("Arial Bold.ttf" if bold else "Arial.ttf"),
                Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold
                     else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"),
                Path(os.environ.get("WINDIR","C:/Windows"))/"Fonts"/("arialbd.ttf" if bold else "arial.ttf")]
    for path in candidates:
        if path.exists():
            return ImageFont.truetype(str(path),size)
    raise RuntimeError("Arial or DejaVu Sans is required.")

def wrap(draw,text,selected_font,width):
    lines=[]
    for paragraph in str(text).splitlines() or [""]:
        current=""
        for word in paragraph.split():
            trial=(current+" "+word).strip()
            if current and draw.textlength(trial,font=selected_font)>width:
                lines.append(current)
                current=word
            else:
                current=trial
        lines.append(current)
    return lines

def capture_path(root,filename):
    path=(Path(root)/filename).resolve()
    if not path.is_relative_to(Path(root).resolve()):
        raise ValueError("Capture escaped its directory.")
    if not path.is_file():
        raise FileNotFoundError("Missing actual-UI capture: "+str(path))
    return path

def compose_frame(capture,keyframe,scene,index,total,caption,speed,width,height,progress):
    canvas=Image.new("RGB",(width,height),C["bg"])
    draw=ImageDraw.Draw(canvas)
    scale=width/1920
    margin,header,footer=round(46*scale),round(120*scale),round(185*scale)
    title_font,small_font,caption_font=font(round(38*scale),True),font(max(12,round(19*scale))),font(round(35*scale))
    draw.text((margin,round(20*scale)),"CHROMATIN STUDIO  /  GUIDED WALKTHROUGH",font=small_font,fill=C["blue"])
    title=f"{index:02d}  {scene['title']}"
    if draw.textlength(title,font=title_font)>width-2*margin-round(180*scale):
        title_font=font(round(32*scale),True)
    draw.text((margin,round(52*scale)),title,font=title_font,fill=C["text"])
    badge="AT YOUR PACE" if speed == 1 else f"{speed:g}×  CAPTION TOUR"
    draw.text((width-margin-draw.textlength(badge,font=small_font),round(28*scale)),badge,font=small_font,fill=C["coral"])
    draw.line((margin,header-round(10*scale),width-margin,header-round(10*scale)),fill=C["line"],width=max(1,round(scale)))
    screenshot=ImageOps.exif_transpose(Image.open(capture)).convert("RGB")
    rect=keyframe.get("focus_rect")
    if rect is not None:
        if not isinstance(rect,list) or len(rect)!=4 or not all(isinstance(v,(int,float)) and math.isfinite(v) for v in rect):
            raise ValueError(scene["id"]+": focus_rect must be [x,y,width,height].")
        x,y,w,h=map(int,rect)
        if x<0 or y<0 or w<=0 or h<=0 or x+w>screenshot.width or y+h>screenshot.height:
            raise ValueError(scene["id"]+": focus rectangle is outside the capture.")
        screenshot=screenshot.crop((x,y,x+w,y+h))
    area=(margin,header,width-margin,height-footer)
    display=ImageOps.contain(screenshot,(area[2]-area[0],area[3]-area[1]),Image.Resampling.LANCZOS)
    left=area[0]+(area[2]-area[0]-display.width)//2
    top=area[1]+(area[3]-area[1]-display.height)//2
    draw.rounded_rectangle((area[0]-1,area[1]-1,area[2]+1,area[3]+1),radius=round(5*scale),fill=C["panel"],outline=C["line"])
    canvas.paste(display,(left,top))
    caption_top=height-footer+round(23*scale)
    lines=wrap(draw,caption,caption_font,width-2*margin-round(60*scale))
    if len(lines)>3:
        caption_font=font(round(30*scale))
        lines=wrap(draw,caption,caption_font,width-2*margin-round(60*scale))
    if len(lines)>3:
        raise ValueError(scene["id"]+": caption is too long for three readable lines.")
    line_height=round(45*scale)
    for row,line in enumerate(lines):
        draw.text((margin+round(16*scale),caption_top+row*line_height),line,font=caption_font,fill=C["text"])
    draw.rectangle((margin,caption_top+round(4*scale),margin+round(3*scale),
                    caption_top+max(round(28*scale),len(lines)*line_height-round(8*scale))),fill=C["coral"])
    py=height-round(22*scale)
    draw.line((margin,py,width-margin,py),fill=C["line"],width=max(2,round(3*scale)))
    draw.line((margin,py,margin+round((width-2*margin)*progress),py),fill=C["blue"],width=max(2,round(3*scale)))
    return canvas

def scene_intervals(scene):
    keyframes,cursor=[],0.
    edges={0.,float(scene["normal_duration_seconds"])}
    for keyframe in scene["keyframes"]:
        end=cursor+float(keyframe["duration"])
        keyframes.append((cursor,end,keyframe))
        edges.update((cursor,end))
        cursor=end
    for caption in scene["captions"]:
        edges.update((float(caption["start"]),float(caption["end"])))
    intervals=[]
    for start,end in zip(sorted(edges),sorted(edges)[1:]):
        if end-start<1e-7:
            continue
        midpoint=(start+end)/2
        keyframe=next(k for a,b,k in keyframes if a<=midpoint<b)
        caption=next((c["text"] for c in scene["captions"] if float(c["start"])<=midpoint<float(c["end"])),"")
        intervals.append((start,end,keyframe,caption))
    return intervals

def concat_quote(path):
    return "'"+str(Path(path).resolve()).replace("'","'\\''")+"'"

def narration_timeline(data, manifest_path, fps, pause):
    """Validate source audio and fit visuals to speech, never speech to visuals."""
    manifest_path=Path(manifest_path).resolve()
    manifest=json.loads(manifest_path.read_text())
    records={entry["id"]:entry for entry in manifest["scenes"]}
    if len(records)!=len(manifest["scenes"]):
        raise ValueError("Narration contains duplicate scene IDs.")
    result=copy.deepcopy(data)
    for scene in result["scenes"]:
        record=records.get(scene["id"])
        if not record or record.get("source_text")!=scene["narration"]:
            raise ValueError(scene["id"]+": narration is missing or stale; regenerate it from this storyboard.")
        path=capture_path(manifest_path.parent,record["file"])
        if sha256(path)!=record["wav_sha256"]:
            raise ValueError(scene["id"]+": narration WAV hash mismatch.")
        with wave.open(str(path),"rb") as audio:
            measured=audio.getnframes()/audio.getframerate()
            expected_bytes=audio.getnframes()*audio.getnchannels()*audio.getsampwidth()
            actual_bytes=0
            while chunk:=audio.readframes(65536):
                actual_bytes+=len(chunk)
            if actual_bytes!=expected_bytes:
                raise ValueError(scene["id"]+": narration WAV payload is truncated.")
        if not math.isfinite(measured) or measured<.5 or not math.isclose(measured,float(record["duration_seconds"]),abs_tol=.001):
            raise ValueError(scene["id"]+": invalid narration duration.")
        if record.get("time_stretch") is not False or record.get("speech_truncated") is not False:
            raise ValueError(scene["id"]+": narration must be unaccelerated and complete.")
        original=float(scene["normal_duration_seconds"])
        duration=math.ceil((measured+pause)*fps)/fps
        factor=duration/original
        for keyframe in scene["keyframes"]:
            keyframe["duration"]*=factor
        for caption in scene["captions"]:
            caption["start"]*=factor
            caption["end"]*=factor
        scene["source_duration_seconds"]=original
        scene["normal_duration_seconds"]=duration
        scene["narration_audio"]={"path":str(path),"sha256":record["wav_sha256"],
            "voice":record["voice"],"engine":manifest["engine"],"duration_seconds":measured,
            "words_per_minute":record["words_per_minute"],"acceleration":1,
            "time_stretch":False,"speech_truncated":False,"reading_pause_seconds":duration-measured}
    result["normal_duration_seconds"]=sum(s["normal_duration_seconds"] for s in result["scenes"])
    result["rendered_duration_seconds"]=result["normal_duration_seconds"]
    result["playback_speed"]=1
    return result

def prepare_audio(scene,work,ffmpeg):
    speech=scene["narration_audio"]
    duration=float(scene["normal_duration_seconds"])
    if duration<speech["duration_seconds"]:
        raise ValueError(scene["id"]+": refusing to truncate speech.")
    output=work/"narration-delivered.wav"
    run([ffmpeg,"-y","-hide_banner","-loglevel","error","-i",speech["path"],"-af","apad",
         "-t",f"{duration:.9f}","-ar","48000","-ac","2","-c:a","pcm_s16le",str(output)])
    return output,speech

def render_scene(scene,index,count,args,ffmpeg,captures,work_root,total_before,total):
    work=work_root/scene["id"];work.mkdir(parents=True,exist_ok=True)
    inputs=[]
    for keyframe in scene["keyframes"]:
        path=capture_path(captures,keyframe["capture"])
        with Image.open(path) as image:
            dimensions=list(image.size)
        inputs.append({"capture":keyframe["capture"],"sha256":sha256(path),"dimensions":dimensions,"focus_rect":keyframe.get("focus_rect")})
    fingerprint=hashlib.sha256(json.dumps({"scene":scene,"inputs":inputs,"index":index,"count":count,"offset":total_before,
        "total":total,"renderer":sha256(__file__),"width":args.width,"height":args.height,"fps":args.fps,
        "silent":args.silent,"speed":args.speed,"crf":args.crf,"preset":args.preset,"threads":args.threads},sort_keys=True).encode()).hexdigest()
    cache_path,clip=work/"scene-manifest.json",work/"scene.mp4"
    if not args.force and cache_path.exists() and clip.exists():
        cached=json.loads(cache_path.read_text())
        if cached.get("fingerprint")==fingerprint:
            print("Reusing "+scene["id"],flush=True)
            return clip,cached
    print("Rendering "+scene["id"],flush=True)
    listing=[]
    for position,(start,end,keyframe,caption) in enumerate(scene_intervals(scene)):
        frame=work/f"frame-{position:03d}.png"
        compose_frame(capture_path(captures,keyframe["capture"]),keyframe,scene,index,count,caption,args.speed,
            args.width,args.height,(total_before+start)/total).save(frame)
        listing.extend(["file "+concat_quote(frame),f"duration {(end-start)/args.speed:.9f}"])
    listing.append("file "+concat_quote(frame))
    sequence=work/"frames.ffconcat";sequence.write_text("ffconcat version 1.0\n"+"\n".join(listing)+"\n")
    audio,speech=(None,{"silent":True}) if args.silent else prepare_audio(scene,work,ffmpeg)
    command=[ffmpeg,"-y","-hide_banner","-loglevel","error","-f","concat","-safe","0","-i",str(sequence)]
    command+=["-i",str(audio),"-map","0:v:0","-map","1:a:0","-c:a","aac","-b:a","160k"] if audio else ["-an"]
    duration=float(scene['normal_duration_seconds'])/args.speed
    command+=["-t",f"{duration:.9f}","-vf",f"fps={args.fps},tpad=stop_mode=clone:stop_duration={duration:.9f}","-r",str(args.fps),"-fps_mode","cfr",
              "-c:v","libx264","-threads",str(args.threads),"-preset",args.preset,"-crf",str(args.crf),"-pix_fmt","yuv420p","-movflags","+faststart",str(clip)]
    run(command)
    manifest={"id":scene["id"],"fingerprint":fingerprint,"normal_seconds":scene["normal_duration_seconds"],
        "delivered_seconds":float(scene["normal_duration_seconds"])/args.speed,"captures":inputs,"speech":speech,
        "clip_sha256":sha256(clip)}
    cache_path.write_text(json.dumps(manifest,indent=2)+"\n")
    return clip,manifest

def timestamp(seconds,comma=False):
    millis=round(seconds*1000)
    hours,millis=divmod(millis,3600000);minutes,millis=divmod(millis,60000);secs,millis=divmod(millis,1000)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}{',' if comma else '.'}{millis:03d}"

def write_accessibility(data,output,stem,speed):
    srt,vtt,chapters,metadata=[],["WEBVTT",""],[],[";FFMETADATA1"]
    offset,number=0.,1
    for scene in data["scenes"]:
        duration=float(scene["normal_duration_seconds"])/speed
        chapters.append({"id":scene["id"],"title":scene["title"],"start":offset,"end":offset+duration,"narration":scene["narration"]})
        title=str(scene["title"]).replace("\\","\\\\").replace("=","\\=").replace(";","\\;").replace("#","\\#").replace("\n"," ")
        metadata+=["[CHAPTER]","TIMEBASE=1/1000",f"START={round(offset*1000)}",f"END={round((offset+duration)*1000)}","title="+title]
        for caption in scene["captions"]:
            start,end=offset+float(caption["start"])/speed,offset+float(caption["end"])/speed
            srt+=[str(number),timestamp(start,True)+" --> "+timestamp(end,True),caption["text"],""]
            vtt+=[str(number),timestamp(start)+" --> "+timestamp(end),caption["text"],""]
            number+=1
        offset+=duration
    (output/(stem+".srt")).write_text("\n".join(srt))
    (output/(stem+".vtt")).write_text("\n".join(vtt))
    (output/(stem+"-chapters.json")).write_text(json.dumps(chapters,indent=2,ensure_ascii=False)+"\n")
    chapter_file=output/(stem+"-chapters.ffmetadata");chapter_file.write_text("\n".join(metadata)+"\n")
    return chapters,chapter_file

def write_player(data,output,stem,chapters,speed,silent):
    buttons="".join('<button data-start="'+str(c["start"])+'"><time>'+timestamp(c["start"])[:8]+'</time><span>'+
                    html.escape(c["title"])+"</span></button>" for c in chapters)
    transcript="".join("<section><h3>"+html.escape(c["title"])+"</h3><p>"+html.escape(c["narration"])+"</p></section>" for c in chapters)
    page="""<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Chromatin Studio — video tutorial</title><style>
:root{color-scheme:dark}*{box-sizing:border-box}body{margin:0;background:#080b0f;color:#edf1f4;font:15px/1.7 system-ui,sans-serif}main{max-width:1280px;margin:auto;padding:35px 24px}header{display:flex;justify-content:space-between;align-items:baseline;gap:20px}h1{font:36px/1.2 Georgia,serif;margin:5px 0 12px}header small{color:#ee7687;white-space:nowrap}.playback{display:flex;align-items:center;gap:12px;margin-top:14px}.playback select{font:inherit;background:#121b25;color:#edf1f4;border:1px solid #58738c;border-radius:4px;padding:5px 10px}.playback span{color:#a7b5c4;font-size:13px}.intro{color:#a7b5c4;margin:0 0 22px;max-width:900px}video{display:block;width:100%;background:#000;border:1px solid #2b3847;border-radius:6px}.layout{display:grid;grid-template-columns:1fr 320px;gap:24px;margin-top:24px}.links{display:flex;flex-wrap:wrap;gap:18px;margin-top:15px}a{color:#79b7de;text-decoration:none;border-bottom:1px solid #314c62}nav{display:grid;gap:5px;align-content:start}nav button{display:flex;gap:12px;background:#121b25;border:1px solid #2b3847;color:#dce6ef;padding:10px 12px;text-align:left;border-radius:4px;cursor:pointer}nav button:hover,nav button.active{border-color:#79b7de;background:#1b2b3b}time{color:#79b7de;font-variant-numeric:tabular-nums;font-size:11px}button:focus-visible,a:focus-visible,summary:focus-visible{outline:2px solid #79b7de;outline-offset:3px}details{border:1px solid #2b3847;border-radius:5px;padding:16px 20px;background:#101821}summary{cursor:pointer;color:#79b7de}h3{font-size:16px;margin:20px 0 7px}section p,.notes{color:#b2bfcc;font-size:13px}.notes{margin:0 0 20px}@media(max-width:800px){main{padding:22px 15px}.layout{grid-template-columns:1fr}header{display:block}h1{font-size:30px}nav{order:-1}}
</style><main><header><h1>Build, inspect, and test a model.</h1><small>At your pace__SILENT__</small></header><p class="intro">A guided walkthrough using actual Chromatin Studio views, with natural-paced neural narration and time to inspect each step. Start at 1×, choose a slower pace below, or jump to a chapter.</p>
<video controls preload="metadata" poster="__STEM__-poster.jpg"><source src="__STEM__.mp4" type="video/mp4"><track kind="captions" src="__STEM__.vtt" srclang="en" label="English key points">Your browser does not support HTML video.</video>
<div class="playback"><label for="playback-speed">Playback speed</label><select id="playback-speed"><option value="0.75">0.75× · slower</option><option value="0.85">0.85× · relaxed</option><option value="1" selected>1× · natural</option><option value="1.15">1.15×</option><option value="1.25">1.25×</option><option value="1.5">1.5×</option></select><span>Try 0.85× for a little more time to follow along.</span></div><div class="links"><a href="__STEM__.mp4" download>Download video</a><a href="__STEM__.srt" download>Timed key points (SRT)</a><a href="transcript.txt">Transcript</a><a href="storyboard.json">Editable storyboard</a></div><div class="layout"><div><p class="notes">On-screen notes summarize each step; the optional text track contains the same key points. The complete narration is in the transcript. The voice is synthesized locally with Kokoro. Synthetic calibration is teaching data; benchmark agreement is not full biological validation. Exporting a bundle does not demonstrate completed GPU execution.</p><details><summary>Read the full narration transcript</summary>__TRANSCRIPT__</details></div><nav aria-label="Tutorial chapters">__BUTTONS__</nav></div></main><script src="player.js" defer></script></html>"""
    page=page.replace("__SPEED__",f"{speed:g}").replace("__SILENT__"," · captions only" if silent else "").replace("__STEM__",html.escape(stem)).replace("__BUTTONS__",buttons).replace("__TRANSCRIPT__",transcript)
    if silent:
        page=page.replace("with natural-paced neural narration and time to inspect each step", "with on-screen key points and time to inspect each step").replace("The voice is synthesized locally with Kokoro. ", "This export has no audio. ")
    (output/"index.html").write_text(page)
    (output/"player.js").write_text("""'use strict';
const video = document.querySelector('video');
const buttons = [...document.querySelectorAll('nav button')];
const speed = document.querySelector('#playback-speed');
speed.addEventListener('change', () => { video.playbackRate = Number(speed.value); });
video.addEventListener('ratechange', () => {
  const known = [...speed.options].find(option => Number(option.value) === video.playbackRate);
  if (known) speed.value = known.value;
});
function updateChapter() {
  let active = buttons[0];
  for (const button of buttons) {
    if (Number(button.dataset.start) <= video.currentTime) active = button;
  }
  for (const button of buttons) {
    button.classList.toggle('active', button === active);
    if (button === active) button.setAttribute('aria-current', 'true');
    else button.removeAttribute('aria-current');
  }
}
for (const button of buttons) {
  button.addEventListener('click', () => {
    video.currentTime = Number(button.dataset.start);
    video.play().catch(() => {});
    video.setAttribute('tabindex', '-1');
    video.focus({preventScroll: true});
    video.scrollIntoView({behavior: window.matchMedia('(prefers-reduced-motion: reduce)').matches ? 'auto' : 'smooth', block: 'start'});
    updateChapter();
  });
}
video.addEventListener('timeupdate', updateChapter);
video.addEventListener('loadedmetadata', updateChapter);
updateChapter();
""")

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--storyboard",type=Path,default=HERE/"storyboard.json")
    parser.add_argument("--captures",type=Path,default=HERE/"captures")
    parser.add_argument("--output-dir",type=Path,default=HERE)
    parser.add_argument("--work-dir",type=Path)
    parser.add_argument("--speed",type=float,default=1,help="Caption-only speed; narrated videos always use 1.")
    parser.add_argument("--silent",action="store_true")
    parser.add_argument("--narration-manifest",type=Path,default=HERE/"narration-neural/narration-manifest.json")
    parser.add_argument("--reading-pause",type=float,default=1.25,help="Extra seconds after each spoken chapter.")
    parser.add_argument("--scene",action="append",help="Render a scene only. Repeat for multiple scenes.")
    parser.add_argument("--width",type=int,default=1920)
    parser.add_argument("--height",type=int,default=1080)
    parser.add_argument("--fps",type=int,default=30)
    parser.add_argument("--crf",type=int,default=20)
    parser.add_argument("--preset",default="veryfast")
    parser.add_argument("--threads",type=int,default=4)
    parser.add_argument("--force",action="store_true")
    parser.add_argument("--validate-only",action="store_true")
    parser.add_argument("--capture-provenance",type=Path)
    args=parser.parse_args()
    if not 1<=args.threads<=16:
        parser.error("Use 1–16 rendering threads.")
    if not args.silent and args.speed!=1:
        parser.error("Neural narration is delivered at its recorded pace. Use the player's speed control while watching.")
    if not math.isfinite(args.reading_pause) or not .5<=args.reading_pause<=10:
        parser.error("Reading pause must be between 0.5 and 10 seconds.")
    if not math.isfinite(args.speed) or not .25<=args.speed<=4 or args.width<1280 or args.height<720 or args.width%2 or args.height%2 or args.fps not in (24,25,30,60) or not 0<=args.crf<=40:
        parser.error("Use speed .25–4, even dimensions at least 1280×720, fps 24/25/30/60 and CRF 0–40.")
    data=load_storyboard(args.storyboard)
    if not args.silent:
        data=narration_timeline(data,args.narration_manifest,args.fps,args.reading_pause)
    ids={s["id"] for s in data["scenes"]}
    if args.scene and set(args.scene)-ids:
        parser.error("Unknown scene ID: "+", ".join(sorted(set(args.scene)-ids)))
    selected=[s for s in data["scenes"] if not args.scene or s["id"] in args.scene]
    for scene in selected:
        for keyframe in scene["keyframes"]:
            capture_path(args.captures,keyframe["capture"])
    if args.validate_only:
        print(json.dumps({"valid":True,"scenes":len(selected),"captures":sum(len(s["keyframes"]) for s in selected),
                          "delivered_seconds":sum(float(s["normal_duration_seconds"])/args.speed for s in selected)}))
        return 0
    ffmpeg=imageio_ffmpeg.get_ffmpeg_exe()
    tag=hashlib.sha256(str(args.storyboard.resolve()).encode()).hexdigest()[:10]
    work=args.work_dir or Path(tempfile.gettempdir())/("chromatin-tutorial-"+tag)
    work.mkdir(parents=True,exist_ok=True);args.output_dir.mkdir(parents=True,exist_ok=True)
    total=sum(float(s["normal_duration_seconds"]) for s in data["scenes"])
    clips,manifests,offset=[],[],0.
    for index,scene in enumerate(data["scenes"],1):
        if scene in selected:
            clip,manifest=render_scene(scene,index,len(data["scenes"]),args,ffmpeg,args.captures,work,offset,total)
            clips.append(clip);manifests.append(manifest)
        offset+=float(scene["normal_duration_seconds"])
    if args.scene:
        for clip,manifest in zip(clips,manifests):
            destination=args.output_dir/(manifest["id"]+"-"+f"{args.speed:g}x"+".mp4")
            shutil.copy2(clip,destination);print("Scene output: "+str(destination),flush=True)
        print("Scene cache updated. Run without --scene to assemble the tutorial.")
        return 0
    stem="chromatin-studio-tutorial" if not args.silent else "chromatin-studio-tutorial-silent-"+f"{args.speed:g}x"
    chapters,metadata=write_accessibility(data,args.output_dir,stem,args.speed)
    playlist=work/"scenes.ffconcat"
    playlist.write_text("ffconcat version 1.0\n"+"\n".join(
        "file "+concat_quote(p)+"\nduration "+f"{m['delivered_seconds']:.9f}"
        for p,m in zip(clips,manifests))+"\n")
    destination=args.output_dir/(stem+".mp4")
    run([ffmpeg,"-y","-hide_banner","-loglevel","error","-f","concat","-safe","0","-i",str(playlist),
         "-f","ffmetadata","-i",str(metadata),"-map","0","-map_metadata","1","-map_chapters","1",
         "-t",f"{total/args.speed:.9f}","-vf",f"setpts=PTS-STARTPTS,fps={args.fps},tpad=stop_mode=clone:stop_duration={total/args.speed:.9f}",
         "-c:v","libx264","-threads",str(args.threads),"-preset",args.preset,"-crf",str(args.crf),"-pix_fmt","yuv420p","-c:a","copy",
         "-movflags","+faststart",str(destination)])
    first=data["scenes"][0]
    poster=compose_frame(capture_path(args.captures,first["keyframes"][0]["capture"]),first["keyframes"][0],first,1,
        len(data["scenes"]),first["captions"][0]["text"],args.speed,args.width,args.height,0)
    poster.save(args.output_dir/(stem+"-poster.jpg"),quality=92)
    write_player(data,args.output_dir,stem,chapters,args.speed,args.silent)
    record={"schema_version":2,"title":data["title"],"storyboard_sha256":sha256(args.storyboard),
        "renderer_sha256":sha256(__file__),"speed":args.speed,"normal_seconds":total,"delivered_seconds":total/args.speed,
        "resolution":[args.width,args.height],"fps":args.fps,"silent":args.silent,
        "video_sha256":sha256(destination),"scenes":manifests,"capture_provenance":None,
        "narration_manifest_sha256":None if args.silent else sha256(args.narration_manifest),
        "timing":"Visuals follow measured narration plus reading pauses; no speech acceleration." if not args.silent else "Caption-only timeline."}
    if args.capture_provenance:
        record["capture_provenance"]={"filename":args.capture_provenance.name,"sha256":sha256(args.capture_provenance)}
    (args.output_dir/"render-manifest.json").write_text(json.dumps(record,indent=2)+"\n")
    (args.output_dir/"render-timeline.json").write_text(json.dumps(data,indent=2,ensure_ascii=False)+"\n")
    print("Rendered "+str(destination),flush=True)
    print(f"Duration: {total/args.speed:.3f} s; narration retains its recorded pace.",flush=True)
    return 0

if __name__=="__main__":
    try:
        raise SystemExit(main())
    except (ValueError,RuntimeError,FileNotFoundError) as error:
        print("Render stopped: "+str(error),file=sys.stderr)
        raise SystemExit(1)

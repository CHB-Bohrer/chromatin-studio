"""Tutorial assets stay scoped and media requests support chapter seeking."""
import tempfile
import threading
import unittest
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer
from pathlib import Path
from unittest.mock import patch

from chromatin_studio import server


class TutorialServerTests(unittest.TestCase):
    video_stems = ("chromatin-studio-tutorial", "chromatin-studio-tutorial-2x")

    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temporary.name)
        cls.assets = cls.root / "tutorial"
        cls.assets.mkdir()
        cls.video = bytes(range(256)) * 2049  # Spans more than two stream chunks.
        for stem in cls.video_stems:
            (cls.assets / f"{stem}.mp4").write_bytes(cls.video)
            (cls.assets / f"{stem}.vtt").write_text("WEBVTT\n")
            (cls.assets / f"{stem}.srt").write_text("1\n00:00:00,000 --> 00:00:01,000\nWelcome.\n")
            (cls.assets / f"{stem}-poster.jpg").write_bytes(b"\xff\xd8\xff\xd9")
            (cls.assets / f"{stem}-chapters.json").write_text('[{"title":"Welcome","start":0}]')
        (cls.assets / "index.html").write_text('<script src="player.js"></script>')
        (cls.assets / "player.js").write_text("document.title = 'Tutorial';")
        (cls.assets / "private.json").write_text('{"private":true}')
        (cls.root / "outside.txt").write_text("outside")
        (cls.assets / "transcript.txt").symlink_to(cls.root / "outside.txt")
        cls.patch = patch.object(server, "TUTORIAL", cls.assets)
        cls.patch.start()
        cls.http = ThreadingHTTPServer(("127.0.0.1", 0), server.Handler)
        cls.thread = threading.Thread(target=cls.http.serve_forever, daemon=True)
        cls.thread.start()
        cls.base = f"http://127.0.0.1:{cls.http.server_port}"

    @classmethod
    def tearDownClass(cls):
        cls.http.shutdown()
        cls.http.server_close()
        cls.thread.join()
        cls.patch.stop()
        cls.temporary.cleanup()

    def get(self, name, byte_range=None):
        headers = {} if byte_range is None else {"Range": byte_range}
        request = urllib.request.Request(self.base + "/tutorial/" + name, headers=headers)
        try:
            response = urllib.request.urlopen(request, timeout=10)
        except urllib.error.HTTPError as exc:
            response = exc
        with response:
            return response.status, response.headers, response.read()

    def test_page_script_captions_and_full_video(self):
        for name in ("index.html", "player.js", "chromatin-studio-tutorial.vtt"):
            status, headers, content = self.get(name)
            self.assertEqual(status, 200)
            self.assertEqual(int(headers["Content-Length"]), len(content))
            self.assertIn("script-src 'self'", headers["Content-Security-Policy"])
            self.assertEqual(headers["X-Content-Type-Options"], "nosniff")
        for stem in self.video_stems:
            with self.subTest(video=stem):
                status, headers, content = self.get(f"{stem}.mp4")
                self.assertEqual(status, 200)
                self.assertEqual(content, self.video)
                self.assertEqual(headers.get_content_type(), "video/mp4")
                self.assertEqual(headers["Accept-Ranges"], "bytes")

    def test_narrated_and_legacy_poster_captions_and_chapters_are_served(self):
        for stem in self.video_stems:
            for suffix, content_type in ((".vtt", "text/vtt"), (".srt", "application/x-subrip"),
                                         ("-poster.jpg", "image/jpeg"),
                                         ("-chapters.json", "application/json")):
                name = stem + suffix
                with self.subTest(asset=name):
                    status, headers, content = self.get(name)
                    self.assertEqual(status, 200)
                    self.assertEqual(headers.get_content_type(), content_type)
                    self.assertEqual(content, (self.assets / name).read_bytes())
                    self.assertEqual(int(headers["Content-Length"]), len(content))
                    self.assertEqual(headers["X-Content-Type-Options"], "nosniff")
                    self.assertIn("script-src 'self'", headers["Content-Security-Policy"])

    def test_video_bounded_open_ended_suffix_and_clamped_ranges(self):
        size = len(self.video)
        for requested, start, end in [("bytes=0-1", 0, 1), ("bytes=262140-262150", 262140, 262150),
                                      ("bytes=500000-", 500000, size-1), ("bytes=-10", size-10, size-1),
                                      (f"bytes={size-10}-{size+999}", size-10, size-1),
                                      (f"bytes=-{size+1}", 0, size-1)]:
            for stem in self.video_stems:
                with self.subTest(video=stem, range=requested):
                    status, headers, content = self.get(f"{stem}.mp4", requested)
                    self.assertEqual(status, 206)
                    self.assertEqual(content, self.video[start:end+1])
                    self.assertEqual(headers["Content-Range"], f"bytes {start}-{end}/{size}")
                    self.assertEqual(int(headers["Content-Length"]), len(content))

    def test_invalid_and_unsatisfiable_ranges_return_416(self):
        for requested in ("bytes=", "bytes=-", "bytes=-0", "bytes=8-2", "bytes=a-b", "items=0-1",
                          "bytes=0-1,4-5", f"bytes={len(self.video)}-", "bytes=" + "9"*100 + "-"):
            for stem in self.video_stems:
                with self.subTest(video=stem, range=requested):
                    status, headers, content = self.get(f"{stem}.mp4", requested)
                    self.assertEqual(status, 416)
                    self.assertEqual(headers["Content-Range"], f"bytes */{len(self.video)}")
                    self.assertEqual(content, b"")

    def test_unpublished_assets_traversal_and_outside_symlink_are_refused(self):
        for name in ("private.json", "../outside.txt", "%2e%2e/outside.txt", "%252e%252e/outside.txt",
                     "subdir/index.html", "index.html%00", "transcript.txt", "render_tutorial.py"):
            with self.subTest(path=name):
                self.assertEqual(self.get(name)[0], 404)


if __name__ == "__main__":
    unittest.main()

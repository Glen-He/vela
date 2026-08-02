from __future__ import annotations

import threading
from collections.abc import Generator
from contextlib import contextmanager
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import ClassVar

from vela.preparation.receptors.acquisition import download_file
from vela.preparation.receptors.models import DownloadConfig, ReceptorError


class _RetryHandler(BaseHTTPRequestHandler):
    attempts: ClassVar[int] = 0
    payload: ClassVar[bytes] = b"validated-payload"

    def do_GET(self) -> None:
        type(self).attempts += 1
        if type(self).attempts == 1:
            self.send_response(503)
            self.send_header("Retry-After", "0")
            self.end_headers()
            return
        self.send_response(200)
        self.send_header("Content-Length", str(len(self.payload)))
        self.send_header("ETag", '"test-etag"')
        self.end_headers()
        self.wfile.write(self.payload)

    def log_message(self, format: str, *args: object) -> None:
        del format, args


@contextmanager
def _retry_server() -> Generator[str]:
    _RetryHandler.attempts = 0
    server = ThreadingHTTPServer(("127.0.0.1", 0), _RetryHandler)
    address = server.server_address
    host = str(address[0])
    port = int(address[1])
    thread = threading.Thread(target=server.serve_forever)
    thread.start()
    try:
        yield f"http://{host}:{port}/structure.cif"
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


def _settings() -> DownloadConfig:
    return DownloadConfig(
        coordinate_base_url="https://coordinates.example.test",
        metadata_base_url="https://metadata.example.test",
        retries=3,
        timeout_seconds=5.0,
        backoff_initial_seconds=0.0,
        backoff_multiplier=2.0,
        chunk_size_bytes=4,
        user_agent="vela-test",
    )


def _validate(path: Path) -> None:
    if path.read_bytes() != _RetryHandler.payload:
        raise ReceptorError("unexpected test payload")


def test_download_retries_then_atomically_installs_and_reuses_cache(
    tmp_path: Path,
) -> None:
    destination = tmp_path / "raw" / "test.cif"
    with _retry_server() as url:
        first = download_file(
            url=url,
            destination=destination,
            settings=_settings(),
            validator=_validate,
        )
        assert first.status == "downloaded"
        assert _RetryHandler.attempts == 2
        assert destination.read_bytes() == _RetryHandler.payload
        assert list(destination.parent.glob("*.part")) == []
        assert list(destination.parent.glob(".*.part")) == []

        second = download_file(
            url=url,
            destination=destination,
            settings=_settings(),
            validator=_validate,
        )
        assert second.status == "cached"
        assert _RetryHandler.attempts == 2

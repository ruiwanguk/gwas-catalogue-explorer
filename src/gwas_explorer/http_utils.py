"""Shared HTTP utilities: retry session and rate-limited executor."""

import logging
import threading
import time
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

logger = logging.getLogger(__name__)

_RETRY_STATUS_CODES = [429, 500, 502, 503, 504]


def create_session(retries: int = 3, backoff_factor: float = 1.0) -> requests.Session:
    """Create a requests session with automatic retry and exponential backoff."""
    session = requests.Session()
    retry = Retry(
        total=retries,
        backoff_factor=backoff_factor,
        status_forcelist=_RETRY_STATUS_CODES,
        allowed_methods=["GET", "POST"],
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    return session


class RateLimitedExecutor:
    """ThreadPoolExecutor that enforces a minimum delay between task starts."""

    def __init__(self, max_workers: int = 4, min_delay: float = 0.1):
        self.max_workers = max_workers
        self.min_delay = min_delay
        self._lock = threading.Lock()
        self._last_start = 0.0

    def _throttled_call(self, fn: Callable, *args):
        with self._lock:
            now = time.monotonic()
            wait = self.min_delay - (now - self._last_start)
            if wait > 0:
                time.sleep(wait)
            self._last_start = time.monotonic()
        return fn(*args)

    def map(self, fn: Callable, items: list[tuple]) -> list:
        """Run fn(*item) for each item, returning results in arbitrary order."""
        results = []
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {
                executor.submit(self._throttled_call, fn, *item): i
                for i, item in enumerate(items)
            }
            for future in as_completed(futures):
                results.append(future.result())
        return results

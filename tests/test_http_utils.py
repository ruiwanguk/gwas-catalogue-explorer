"""Tests for shared HTTP utilities: retry session and rate-limited executor."""

import time

import requests
import responses

from gwas_explorer.http_utils import RateLimitedExecutor, create_session


class TestCreateSession:
    """Tests for create_session factory."""

    def test_returns_session(self) -> None:
        session = create_session()
        assert isinstance(session, requests.Session)

    def test_mounts_retry_adapters(self) -> None:
        session = create_session(retries=5, backoff_factor=2.0)
        for prefix in ("http://", "https://"):
            adapter = session.get_adapter(prefix)
            assert adapter.max_retries.total == 5
            assert adapter.max_retries.backoff_factor == 2.0

    def test_retry_status_codes(self) -> None:
        session = create_session()
        adapter = session.get_adapter("https://")
        assert set(adapter.max_retries.status_forcelist) == {429, 500, 502, 503, 504}

    def test_allowed_methods(self) -> None:
        session = create_session()
        adapter = session.get_adapter("https://")
        assert "GET" in adapter.max_retries.allowed_methods
        assert "POST" in adapter.max_retries.allowed_methods

    @responses.activate
    def test_successful_request(self) -> None:
        responses.add(responses.GET, "https://example.com/api", json={"ok": True}, status=200)
        session = create_session()
        resp = session.get("https://example.com/api")
        assert resp.status_code == 200
        assert resp.json() == {"ok": True}


class TestRateLimitedExecutor:
    """Tests for RateLimitedExecutor."""

    def test_map_returns_all_results(self) -> None:
        executor = RateLimitedExecutor(max_workers=2, min_delay=0.0)
        results = executor.map(lambda x, y: x + y, [(1, 2), (3, 4), (5, 6)])
        assert sorted(results) == [3, 7, 11]

    def test_map_empty_items(self) -> None:
        executor = RateLimitedExecutor(max_workers=2, min_delay=0.0)
        results = executor.map(lambda x: x, [])
        assert results == []

    def test_map_single_item(self) -> None:
        executor = RateLimitedExecutor(max_workers=1, min_delay=0.0)
        results = executor.map(lambda x: x * 2, [(5,)])
        assert results == [10]

    def test_min_delay_enforced(self) -> None:
        """Tasks should be spaced by at least min_delay."""
        timestamps: list[float] = []

        def record_time(x: int) -> int:
            timestamps.append(time.monotonic())
            return x

        executor = RateLimitedExecutor(max_workers=1, min_delay=0.05)
        executor.map(record_time, [(1,), (2,), (3,)])

        timestamps.sort()
        for i in range(1, len(timestamps)):
            gap = timestamps[i] - timestamps[i - 1]
            assert gap >= 0.04, f"Gap {gap:.4f}s is less than min_delay"

    def test_propagates_exceptions(self) -> None:
        def fail(x: int) -> int:
            raise ValueError("boom")

        executor = RateLimitedExecutor(max_workers=1, min_delay=0.0)
        try:
            executor.map(fail, [(1,)])
            assert False, "Should have raised"
        except ValueError as e:
            assert "boom" in str(e)

from __future__ import annotations

from molpy.core.buckets import BucketRegistry, TypeBucket


class TestTypeBucket:
    def test_add_discard_iter(self) -> None:
        b: TypeBucket[int] = TypeBucket()
        b.add(1, 2)
        assert set(iter(b)) == {1, 2}
        b.discard(2)
        assert set(iter(b)) == {1}


class TestBucketRegistry:
    class A: ...
    class B(A): ...

    def test_register_and_add(self) -> None:
        reg = BucketRegistry()
        reg.register(self.A)
        reg.register(self.B)
        a = self.A()
        b = self.B()
        reg.add_by_instance(a)
        reg.add_by_instance(b)
        assert a in set(reg.bucket(self.A))
        assert b in set(reg.bucket(self.B))
        reg.discard_any(a)
        assert a not in set(reg.bucket(self.A))
        assert self.A in reg.classes() and self.B in reg.classes()

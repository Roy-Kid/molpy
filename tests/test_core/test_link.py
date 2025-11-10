from molpy.core.entity import Entity, Link

class TestLink:
    def test_endpoints_and_replace(self) -> None:
        a = Entity({"id": 1})
        b = Entity({"id": 2})
        c = Entity({"id": 3})
        l = Link([a, b], order=1)
        assert l.endpoints == (a, b)
        assert l["order"] == 1
        l.replace_endpoint(b, c)
        assert l.endpoints == (a, c)

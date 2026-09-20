"""NotPublic: a class attribute that names the public door instead."""

import pytest

from molpy.core.entity import NotPublic


class _Graph:
    add_node = NotPublic("def_node")

    def def_node(self) -> str:
        return "public"


def test_reading_the_hidden_name_names_the_door():
    with pytest.raises(AttributeError, match="add_node is not public; use def_node"):
        _Graph().add_node


def test_only_the_hidden_name_is_affected():
    assert _Graph().def_node() == "public"
    assert not hasattr(_Graph(), "add_node")

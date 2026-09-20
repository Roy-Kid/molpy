"""Capabilities molpy publishes for other molcrafts products to consume.

Nothing here is imported by molpy itself, and nothing here imports another
molcrafts product. Each module implements a structural contract that a host
looks up by entry point, so the dependency arrow points one way: molpy owns
the formats and publishes what it can do with them; a platform or a viewer
consumes that without molpy knowing either exists.
"""

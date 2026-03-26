# Developer Guide

This guide is for contributors and library developers. It covers conventions, extension patterns, and internal architecture.

## Part 1: Conventions

How the project works day-to-day.

- [Development Setup](development-setup.md) — clone, install, run tests
- [Contributing](contributing.md) — workflow, PR checklist
- [Coding Style](coding-style.md) — naming, formatting, immutability rules
- [Testing](testing.md) — pytest conventions, markers, coverage
- [Release Process](release-process.md) — versioning, tagging, CI

## Part 2: Extending the Tool Layer

How to add new capabilities using MolPy's established extension points. These patterns have stable interfaces and clear contracts — you implement a subclass or register a handler.

- [Adding a Tool or Compute Operation](extending-tools.md) — `Tool`, `Compute`, custom recipes
- [Adding an I/O Format](extending-io.md) — readers, writers, force field formatters
- [Adding a Wrapper or Adapter](extending-integration.md) — external CLI tools, in-memory bridges

## Part 3: Extending Core Structures

How to extend MolPy's foundational data model. These changes are deeper and require understanding the type bucket system, the force field hierarchy, and the potential dispatch chain.

- [Extending the Data Model](extending-core.md) — new Entity/Link types, custom Struct subclasses
- [Extending the Force Field](extending-forcefield.md) — custom Style, Type, Potential, and formatter registration

## Integrations

- [MCP for LLM Agents](../user-guide/mcp.md) — exposing MolPy's source to LLM agents via Model Context Protocol


## Support

- Issues: <https://github.com/MolCrafts/molpy/issues>
- Discussions: <https://github.com/MolCrafts/molpy/discussions>

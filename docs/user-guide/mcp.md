# MCP: Letting LLM Agents Read MolPy's Source

## Agents need to understand the library, not just call it

An LLM agent working with MolPy faces a fundamental problem: its training data may be outdated or incomplete. It can guess at API names and parameter types, but those guesses break whenever the library changes. Even when the guess is correct, the agent has no way to verify it — it cannot read the source code to confirm that a function exists, what it accepts, or what it returns.

**MolPy's MCP server solves this by exposing the library's source code, docstrings, and signatures as structured tool calls over the Model Context Protocol.** Instead of guessing, the agent queries the server: "What symbols does `molpy.core.atomistic` export?" or "Show me the source of `Reacter.run`." The server returns accurate, up-to-date answers because it reads directly from the installed package.

This is not a remote execution server. The agent does not call `PrepareMonomer` through MCP. Instead, it uses MCP to *understand* MolPy's API, then writes Python code that calls the API directly. The MCP server is a documentation retrieval tool, not a computation engine.

## Six tools for code exploration

The server exposes six tools. Together they let an agent navigate the package tree, understand any symbol, and search the codebase.

**`list_modules`** returns all importable module paths under a given prefix. An agent starting from scratch calls `list_modules("molpy")` to discover the package structure — `molpy.core.atomistic`, `molpy.tool.polymer`, `molpy.parser.smiles`, and so on.

**`list_symbols`** returns the public symbols exported by a module, each with a one-line summary from its docstring. This tells the agent what a module contains without reading the full source.

**`get_docstring`** retrieves the cleaned docstring of any module, class, or function. Since MolPy docstrings follow Google style with `Args`, `Returns`, `Preferred for`, and `Avoid when` sections, the agent gets structured usage guidance directly from the source.

**`get_signature`** returns the call signature of a callable — parameter names, types, and defaults. Combined with the docstring, this gives the agent everything it needs to write a correct function call.

**`get_source`** retrieves full source code when the agent needs to understand implementation details — for example, what steps `PrepareMonomer.run()` actually performs, or how `Connector` detects leaving groups.

**`search_source`** does a case-insensitive substring search across all `.py` files, returning file paths, line numbers, and matching lines. This is the agent's grep — useful for finding where a class is defined, where a function is called, or how an error message originates.

## Installation and setup

MCP support requires the `fastmcp` dependency:

```bash
pip install molpy[mcp]
```

Start the server from the CLI:

```bash
molpy mcp                            # stdio transport (default)
molpy mcp -t streamable-http         # HTTP on 127.0.0.1:8787
molpy mcp -t streamable-http -p 9000 # custom port
```

To connect from Claude Code, add the server to your MCP config:

```json
{
  "mcpServers": {
    "molpy": {
      "command": "molpy",
      "args": ["mcp"]
    }
  }
}
```

For Cursor or other HTTP-based clients, start the server with `streamable-http` transport and point the client at `http://127.0.0.1:8787`.

## How an agent uses the server in practice

A typical interaction looks like this. The agent receives a task — "build a 10-unit PEO chain and export to LAMMPS." It does not know the current API, so it starts by exploring:

1. `list_modules("molpy.tool")` — discovers `molpy.tool.polymer`.
2. `list_symbols("molpy.tool.polymer")` — sees `PrepareMonomer`, `polymer`, `BuildPolymer`, etc.
3. `get_docstring("molpy.tool.polymer.polymer")` — reads the usage guidance, learns about notation auto-detection.
4. `get_signature("molpy.tool.polymer.polymer")` — gets the exact parameter names and defaults.

With this information, the agent writes correct Python code:

```python
from molpy.tool import polymer
chain = polymer("{[<]CCO[>]}|10|")
```

The key difference from guessing: every API call the agent writes is grounded in the actual source. If MolPy renames a parameter or moves a function, the agent discovers the change through MCP rather than failing silently with stale knowledge.

## See Also

- [Tool Layer](tools.md) — what the Tool recipes do and when to use them
- [API Reference: Tool](../api/tool.md) — full parameter documentation

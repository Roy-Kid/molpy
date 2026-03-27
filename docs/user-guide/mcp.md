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


## Installation and adding the server

MCP support requires the `fastmcp` dependency:

```bash
pip install molpy[mcp]
```

### Claude Code

One command registers the server. The `--scope` flag controls whether it applies to the current project or to all your projects:

```bash
# Project-level (recommended) — adds to .mcp.json in the repo root
claude mcp add molpy --scope project -- molpy mcp

# User-level — adds to ~/.claude/settings.json
claude mcp add molpy -- molpy mcp
```

That's it. Start a new Claude Code session and the server is available. To verify, type `/mcp` in the prompt — you should see `molpy` with six tools listed.

!!! tip "Virtual environments"
    If MolPy is installed inside a virtual environment and `molpy` is not on the system PATH, pass the full path:

    ```bash
    claude mcp add molpy -- /path/to/venv/bin/molpy mcp
    ```

    Or use `uv run` to handle activation automatically:

    ```bash
    claude mcp add molpy -- uv run --directory /path/to/molpy molpy mcp
    ```

### Claude Desktop

Claude Desktop reads from its own config file:

- **macOS**: `~/Library/Application Support/Claude/claude_desktop_config.json`
- **Windows**: `%APPDATA%\Claude\claude_desktop_config.json`

Add a `mcpServers` entry:

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

Restart Claude Desktop after saving. The MolPy tools appear in the tool picker (the hammer icon).

## How MCP works

MCP (Model Context Protocol) is a lightweight RPC layer between an LLM client and a tool server. Understanding the protocol helps you debug connection issues and reason about what the agent is doing.

### The client-server model

The MolPy MCP server is a **separate process** that the LLM client launches and communicates with. The data flow looks like this:

```text
┌─────────────┐   stdin/stdout    ┌──────────────────┐
│  LLM Client │ ◄──────────────► │  molpy mcp        │
│  (Claude)   │   JSON-RPC msgs   │  (MCP Server)     │
└─────────────┘                   └──────────────────┘
       │                                   │
       │  "What tools do you have?"        │
       │ ─────────────────────────────────►│
       │                                   │
       │  "list_modules, list_symbols, .." │
       │ ◄─────────────────────────────────│
       │                                   │
       │  call list_modules("molpy.tool")  │
       │ ─────────────────────────────────►│
       │                                   │
       │  ["molpy.tool", "molpy.tool.base",│
       │   "molpy.tool.polymer", ...]      │
       │ ◄─────────────────────────────────│
```
When the client starts, it sends a handshake request asking the server to list its capabilities. The server responds with six tool definitions — names, parameter schemas, and descriptions. The LLM sees these definitions as available tools, just like file reading or web search. When the LLM decides to call a tool, the client serializes the call as a JSON-RPC message, sends it to the server over stdin, and reads the result from stdout.

The `stdio` transport is the default and works everywhere. The `streamable-http` and `sse` transports serve the same protocol over HTTP, which is useful for clients that cannot launch subprocesses. The primitive way to call mcp server via agent looks like this:

```python
"Please list all the modules included in molpy."
  ┌─────────────────┬─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
  │     Package     │                                                         Description                                                         │
  ├─────────────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ molpy.adapter   │ Adapters for OpenBabel, RDKit                                                                                               │
  ├─────────────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ molpy.builder   │ Crystal and polymer builders (AmberTools integration, stochastic generation, sequences)                                     │
  ├─────────────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ molpy.cli       │ Command-line interface                                                                                                      │
  ├─────────────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ molpy.compute   │ Computations (MCD, PMSD, RDKit-based, time series)                                                                          │
  ├─────────────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ molpy.core      │ Core data structures — atomistic/CG models, box, element, entity, forcefield, frame, topology, trajectory, units, selectors │
  ├─────────────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ molpy.data      │ Data resources (forcefield data)                                                                                            │
  ├─────────────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ molpy.engine    │ Simulation engines (CP2K, LAMMPS)                                                                                           │
  ├─────────────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ molpy.io        │ I/O for data formats (AC, Amber, GRO, H5, LAMMPS, MOL2, PDB, TOP, XSF, XYZ), forcefields (Amber, frcmod, LAMMPS,            │
  │                 │ Moltemplate, TOP, XML), trajectories (H5, LAMMPS, XYZ), and stores (H5, Zarr)                                               │
  ├─────────────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ molpy.op        │ Operations (geometry)                                                                                                       │
  ├─────────────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ molpy.optimize  │ Optimization (L-BFGS, potential wrappers)                                                                                   │
  ├─────────────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ molpy.pack      │ Packing (constraints, targets, Packmol integration)                                                                         │
  ├─────────────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ molpy.parser    │ Parsers for SMARTS and SMILES (including BigSMILES, cgSMILES, gBigSMILES)                                                   │
  ├─────────────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ molpy.potential │ Potentials — bond, angle, dihedral, improper, pair (LJ, Coulomb)                                                            │
  ├─────────────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ molpy.reacter   │ Reaction tools (connectors, selectors, templates, topology detection)                                                       │
  ├─────────────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ molpy.tool      │ Analysis tools (cross-correlation, MSD, polymer, RDKit, time series)                                                        │
  ├─────────────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ molpy.typifier  │ Force field typing (GAFF, OPLS, graph-based matching, layered engine)                                                       │
  ├─────────────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
  │ molpy.wrapper   │ External tool wrappers (antechamber, prepgen, tleap)                                                                        │
  └─────────────────┴─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘

  Total: 148 modules across 18 top-level packages.

```

Writing effective prompts

The MCP server gives the agent access to MolPy’s API — but it does not guarantee a correct result. The difference between a working script and a broken one almost always comes down to how you phrase the prompt.

A good prompt makes the task fully specified and unambiguous. A bad prompt leaves the agent guessing.

Describe what you want to build — not how to build it

Focus on the system, not the API.

The agent can discover functions through MCP. If you tell it which function to call, you are bypassing that mechanism and increasing the chance of errors.

Weak prompt	Better prompt
Use polymer() to build a PEG chain	Build a PEG chain with 15 repeat units
Call Molpack to pack molecules	Pack 15 chains into a 20 nm cubic box
Use the Box class	Create a periodic simulation box for the system

👉 Rule of thumb:
If your prompt mentions function names, it is probably too low-level.

Always include the key physical parameters

Molecular systems are defined by numbers. If you omit them, the agent has to guess — and guesses are usually wrong.

At minimum, specify:

molecule type (e.g. PEG, PEO, polystyrene)
chain length
number of molecules
box size or density
output format (if needed)

Good example:

Generate a 20 × 20 × 20 nm box containing 15 PEG chains,
each with a length of 15 monomers. Export to LAMMPS DATA format.

This prompt is complete. The agent can start immediately without asking questions.

Keep one prompt = one workflow

Do not combine everything into one request.

❌ Bad:

Build a polymer system, run MD, and compute RDF.

✔ Better:

1. Build 15 PEG chains of length 15 and pack into a 20 nm box
2. Set up a LAMMPS equilibration at 300 K
3. Compute the radial distribution function

👉 Think like a modeling workflow, not a single command.

State constraints that affect the result

If something changes how the system is built, say it explicitly.

Examples:

Use the Amber backend with GAFF2 parameters
Do not use RDKit

If you don’t specify constraints, the agent will pick defaults — which may not match your setup.

Let the agent explore (don’t over-guide it)

The agent will call MCP tools like list_modules or get_signature to understand the API. This is expected.

Avoid:

telling it which function to use
pasting code snippets from memory
forcing a specific implementation

If it makes a mistake after exploring, that usually means:

the API is unclear, or
the documentation needs improvement

That’s useful feedback.

Quick checklist

Before sending a prompt, ask yourself:

Does this describe a real system, not API calls?
Are all numbers specified?
Is this one task, not three?
Did I include constraints that matter?

If yes, the agent will usually get it right on the first try.

## Worked example: 15 PEG chains in a 20 nm box

The rest of this page shows a real session with Claude Code. The MolPy MCP server is connected; the user types one sentence and Claude explores the API, writes a script, and runs it.

### The prompt

```text
Generate a 20 x 20 x 20 nm box containing 15 PEG chains, each of
length 15. Export to LAMMPS DATA format.
```

### Agent exploration

Claude starts by discovering what MolPy offers. The MCP tool calls it makes, and the key information it extracts:

<!-- TODO: Paste actual Claude MCP tool call sequence here.
     Expected flow:
     1. list_modules("molpy") → discovers molpy.tool, molpy.pack, molpy.io, etc.
     2. list_symbols("molpy.tool.polymer") → finds polymer(), polymer_system()
     3. get_signature("molpy.tool.polymer.polymer") → learns the parameter interface
     4. get_docstring("molpy.tool.polymer.polymer") → reads G-BigSMILES notation
     5. list_symbols("molpy.pack") → finds Molpack, InsideBoxConstraint
     6. get_signature("molpy.pack.Molpack.optimize") → learns packing interface
     7. get_signature("molpy.core.box.Box.cubic") → learns box creation
     8. search_source("write_lammps") → finds export function
-->

### The generated script

<!-- TODO: Paste the Python script that Claude produces.
     Expected to cover:
     - import molpy as mp
     - Build 15 chains with mp.tool.polymer("{[<]OCCO[>]}|15|")
     - Pack with Molpack + InsideBoxConstraint(length=200 A)
     - Attach Box.cubic(200.0) to metadata
     - Export with mp.io.write_lammps_data(...)
-->

### Output

<!-- TODO: Paste terminal output from running the script.
     Expected:
     - "built 15 chains, N atoms each"
     - "packed: M total atoms in 200 A cubic box"
     - "wrote peg_box/peg_15x15.data"
     If Packmol is missing, show the error and add a note:
     Packmol must be installed separately — see https://...
-->


## See Also

- [Tool Layer](tools.md) — what the Tool recipes do and when to use them
- [Polydisperse Systems](05_polydisperse_systems.md) — end-to-end workflow from distribution to LAMMPS export
- [API Reference: Tool](../api/tool.md) — full parameter documentation

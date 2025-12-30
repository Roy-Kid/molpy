# LLM-friendly usage

## 1. Why MolPy is LLM-friendly

MolPy assumes that a workflow can be expressed as explicit data plus explicit transformations. Instead of hiding state in a long-lived simulation object, MolPy encourages you to build and modify concrete structures (for example an `Atomistic` graph or a tabular `Frame`) that can be inspected, serialized, and passed between steps.

This matters for LLM-based agents because their “working memory” is limited to what is visible in the prompt and tool results. When an operation’s inputs and outputs are explicit, an agent can reason about what it is changing, verify intermediate results, and repeat a step without relying on unstated context.

MolPy also favors deterministic, inspectable transformations. Many operations are “given X, produce Y” with minimal side effects, so an agent can re-run them and compare outcomes. If a result looks wrong, the agent can inspect the current structure (counts, fields, topology) and decide what to adjust.

Finally, MolPy separates modeling logic from execution backends. The core modeling layer focuses on representing chemistry and structured data. Execution-heavy work (external tools, engines, binaries) is kept in integration layers. This separation lets an agent build and validate a system in-process before asking an external engine to run it.

## 2. Context: MCP support in MolPy

MCP (Model Context Protocol) is a standard way to expose software capabilities as **tool calls** with **machine-readable schemas**. In practice, MCP turns “a Python library with many functions” into “a toolbox” that an LLM agent can call reliably: each tool has a name, a defined input shape, and a structured output.

We use MCP because a plain LLM is good at *planning* but unreliable at *executing* arbitrary code. Without MCP, an agent tends to guess APIs, hallucinate parameters, or have to dive into the source code. MCP provides a strict boundary: the agent chooses *what to do next* and supplies parameters, while MolPy executes the operation and returns structured results that can be inspected and validated. This makes the workflow reproducible, debuggable, and less sensitive to prompt noise.

When MolPy is used by an LLM-based agent, MCP is used to constrain how the agent interacts with the library.
Concrete MolPy operations are exposed as callable tools, while [Context7](https://context7.com/) is used to retrieve the relevant MolPy documentation before any operation is executed. This ensures that tool calls are based on the current documented usage rather than on assumptions or prior knowledge.

To let an LLM agent consult the MolPy documentation before acting, you need to install and enable [Context7](https://context7.com/) in your own editor. [Context7](https://context7.com/) is not part of MolPy and is not installed automatically. Instead, it is configured at the agent or MCP host level, where it exposes documentation retrieval as a callable tool.

Your prompt should enforce three things: (i) always consult [Context7](https://context7.com/) before using MolPy APIs, (ii) copy the exact API names from the excerpt, and (iii) return a minimal plan + the next tool call.

A good prompt pattern is:

1. State the task in one sentence.
2. Require a [Context7](https://context7.com/) lookup first.
3. Require a short plan and then execute via MCP tools.

Here is a concrete template you can reuse:

```text
Task: <one sentence describing what you want to do in MolPy>.

Constraints:
- <constraint 1>

Rules:
- You MUST use Context7 to retrieve the relevant MolPy documentation before writing any code.
- You MUST base your plan and code strictly on the retrieved documentation (do NOT guess APIs).
- If the documentation is insufficient or unclear, you MUST stop and state what is missing.

Steps:
1) Use Context7 to fetch guidance for "<topic>".
2) Provide a short plan that satisfies the Constraints and references.
3) Implement the code following the plan. Keep it minimal and runnable.
```

This prompt style prevents the agent from guessing. It forces an explicit “docs → plan → execute” loop, which is the whole point of adding MCP + Context7 in the first place.

Let's start with an example using MolPy to build a water box with tip3p force field. The prompt would be:

```text
Task: Build a periodic TIP3P water box using MolPy with Packmol, then export LAMMPS force-field files and a LAMMPS data file.

Constraints:
- Use TIP3P water.
- Use Packmol for packing (no manual placement).
- Use a cubic periodic box.
- Export artifacts needed to run LAMMPS: (1) data file and (2) force-field/parameter files required by the chosen MolPy export path.

Rules:
- You MUST use Context7 to retrieve the relevant MolPy documentation before writing any code.
- You MUST base your plan and code strictly on the retrieved documentation (do NOT guess APIs).
- If the documentation is insufficient or unclear, you MUST stop and state what is missing.

Steps:
1) Use Context7 to fetch guidance for "TIP3P water, Packmol packing into a cubic box, converting to Frame, and exporting LAMMPS data + force-field files + input script".
2) Provide a short plan that satisfies the Constraints and references.
3) Implement the code following the plan. Keep it minimal and runnable.
```

Sending this prompt to an LLM-based agent with MCP + Context7 enabled should yield a correct plan and code snippet that accomplishes the task.

## 3. Summary

MolPy’s suitability for LLM-assisted work follows from its design: explicit data models, inspectable transformations, and a clear boundary between modeling and backend execution.

This makes workflows more predictable and easier to reproduce, whether they are driven by a human user, by an LLM-based agent calling tools, or by a combination of both.

MolPy can be used directly as a Python library, and it can also be wrapped behind schema-defined tools (such as via MCP) when you want a stricter, tool-based interaction model.

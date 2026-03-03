# Sample program to check the qasm parser

This is a sample program to check tiny header-only c++ parser for OpenQASM (v2 & a minimal v3 subset). This program
1. reads a QASM file (exported from Python/Qiskit)
2. parses it into a simple IR, and
3. dumps the program structure (registers and instructions).

It's meant for quick experiments and simulation prototypes.

## Requirement

### Build
- C++17 compiler (GCC $\geq$ 9, Clang $\geq$ 10, or MSVC $\geq$ 2019)
- A shell to run the sample build commands

### Generate QASM (Phyton,optional)
- Python 3.9+
- qiskit (for exporting OpenQASM 2.0/3.0)
Install:
```
python3 -m venv path_to_env
source path_to_env/bin/activate
pip3 install qiskit
```

## What's Supported (Quick Summary)
- OpenQASM 2.0
  - `OPENQASM 2.0;`, `include` (ignored), `qreg/creg`,
  - standard gates (`h,x,y,z,s,sdg,t,tdg,rx,ry,rz,u1,u2,u3,cx,cz,swap,ccx,cswap,id`),
  - `barrier`, `reset`, `measure q[i] -> c[i];`,
  - and whole-register expansion (e.g., `h q;`)
- OpenQASM 3.0
  - `OPENQASM 3.0;`, `include` (ignored), `qubit/bit` declarations,
  - assignment-form measurements (`c[i] = measure q[i];`, `cit x = measure q[j];`),
  - standard gate calls (no gate modifiers),
  - `barrier`, `reset`, range slices `q[lo:hi]` (half-open `[lo,hi)`),
  - whole-register expantion.
  - **Note supported**: control flow, gate definitions, timing/pulse, gate modifiers (`ctrl/inv/pow`).

Parsed IR: see `include/qasm/ir.h` (`Program`, `Instruction`, `Op`, etc.)

## Sample Program (C++)
- `sample/qasm_parse/main.cc` (reads QASM from a file or stdin and dumps the IR).

## Build

From repository root (`path_to_home/`):
```
g++ -std=c++17 sample/qasm_parse/main.cc -Iinclude -o parser
# or
clang++ -std=c++17 sample/qasm_parse/main.cc -Iinclude -o parser
```
If you build from `path_to_home/sample/qasm_parse/`:
```
g++ -std=c++17 main.cc -I../../include -o parser
```

## Test 
### Generate a QASM file with Qiskit (OpenQASM 2.0)
Run:
```
cd path_to_home/sample/qasm_parse
python3.13 gen_qasm.py
```

### Run the parser
Pass the file path:
```
# from path_to_home/sample/qasm_parse/
./parser circuit.qasm
```
Or pipe from stdin:
```
# from path_to_home/sample/qasm_parse/
cat circuit.qasm | ./parser
```
You'll see a dump like:
```
qregs:
  qreg q[3]
cregs:
  creg c[3]
instructions: 4
  line 5: 6(h)  q[0]
  line 6: 15(cx)  q[0], q[1]
  line 7: 5(rz) (1.5708)  q[2]
  line 8: 20(measure)  q[1] -> c[1]
```
Note: the numeric op codes correspond to `enum qasm::Op`.
Optionally, add a helper like `op_name(Op)` to print friendly names.

## Notes and Tips
- **Whole-register expansion**: statements line `h q;` expand into per-qubit instructions in IR.
- **Slices (v3)**: `q[lo:hi]` is half-open and expands to `lo, lo+1, ..., hi-1`.
- **Unsupported syntax** in v3 (control flow, gate modifiers `ctrl/inv/pow`, timing/pulse, etc) raises a parse error.
- If you hit `not found` for headers, check your `-I` path points to the parent of `qasm/` (i.e., `../include`), and verify filenames (`.h`) and case sensitivity.

## Extending
- Hook in your simulator by walking `Program.instructions`.

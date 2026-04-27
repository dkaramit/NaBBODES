# NaBBODES Developer Guide

This document is for contributors who want to:

1. understand solver internals,
2. add a new integration **method** in `METHOD.hpp`,
3. hack or customize core algorithm behavior.

It complements `docs/USAGE_AND_API.md` (user-facing API).

---

## 1) High-level architecture

## 1.1 Solver families

- `rkf::Solver<...>` in `RKF/`
- `rosenbrock::Solver<...>` in `Rosenbrock/`

Both are header-only style assemblages:

- class declaration in `*_class.hpp`,
- implementation split across `*_reset.hpp`, `*_calc_k.hpp`, `*_sums.hpp`, step-control files, and `*_steps.hpp`,
- umbrella includes in `RKF/RKF.hpp` and `Rosenbrock/Rosenbrock.hpp`.

## 1.2 Shared design pattern

Both solvers:

- hold method coefficients in template type `RK_method`,
- support compile-time step controller selection via `step_controllers` enum (`simple` / `PI`),
- expose identical parameter struct style (`std::optional` fields),
- store full output history (`time`, `solution`, `error`) and provide accessor API.

---

## 2) Algorithm overview for hackers

## 2.1 RKF (`rkf::Solver`)

Relevant files:

- `RKF/RKF_class.hpp`
- `RKF/RKF_calc_k.hpp`
- `RKF/RKF_sums.hpp`
- `RKF/RKF_step_control_simple.hpp`
- `RKF/RKF_step_control_PI.hpp`
- `RKF/RKF_steps.hpp`

Core flow in `next_step()`:

1. initialize trial step and controller history terms,
2. repeatedly:
   - compute stage values `k` (`calc_k()`),
   - compute solution combination(s) (`sum_bk()` etc),
   - compute embedded error estimate `abs_delta`,
   - apply chosen controller (`step_control()`),
3. once accepted:
   - update state/time,
   - push to output vectors.

Conceptually: this is an adaptive embedded RK integrator with local error control.

## 2.2 Rosenbrock (`rosenbrock::Solver`)

Relevant files:

- `Rosenbrock/Ros_class.hpp`
- `Rosenbrock/Ros_calc_k.hpp`
- `Rosenbrock/Ros_LU.hpp`
- `Rosenbrock/LU/LU.hpp`
- `Rosenbrock/Ros_sums.hpp`
- `Rosenbrock/Ros_step_control_simple.hpp`
- `Rosenbrock/Ros_step_control_PI.hpp`
- `Rosenbrock/Ros_steps.hpp`

Core flow in `next_step()`:

1. evaluate Jacobian + `dfdt` (default finite differences or user callback),
2. build and factor stage linear operator `(I - gamma*h*J)` (LU path),
3. compute stage vectors `k_i`,
4. combine stages into main/embedded updates,
5. estimate local error and run selected step controller,
6. on acceptance, commit state and append outputs.

Conceptually: linearly implicit Rosenbrock method with embedded error estimate and adaptive step size.

---

## 3) Adding a new method in `METHOD.hpp`

This section is the key extension recipe.

## 3.1 Add a new RKF method (`RKF/METHOD.hpp`)

Each RKF method is a struct template:

```cpp
template<class LD>
struct MyMethod {
    static constexpr unsigned int s = ...;   // number of stages
    static constexpr unsigned int p = ...;   // principal order

    using arr  = std::array<LD, s>;
    using arr2 = std::array<std::array<LD, s>, s>;

    static constexpr arr  c = { ... };
    static constexpr arr  b = { ... };
    static constexpr arr  bstar = { ... };   // embedded pair weights
    static constexpr arr2 a = { ... };       // stage matrix
};
```

### Minimum coefficient requirements (RKF)

- `s`, `p`
- `a`, `b`, `bstar`, `c`

### Practical checklist (RKF)

1. Add method struct in `RKF/METHOD.hpp`.
2. Use it in example type alias:
   ```cpp
   using SOLVER = rkf::Solver<n_eqs, MyMethod<LD>, LD, rkf::step_controllers::PI>;
   ```
3. Build:
   ```bash
   make -C RKF/Examples
   ```
4. Run and sanity-check error behavior:
   ```bash
   ./RKF/Examples/Example.run | head
   ```

## 3.2 Add a new Rosenbrock method (`Rosenbrock/METHOD.hpp`)

Each Rosenbrock method is a struct template:

```cpp
template<class LD>
struct MyRosMethod {
    static constexpr unsigned int s = ...;   // stages
    static constexpr unsigned int p = ...;   // order

    using arr  = std::array<LD, s>;
    using arr2 = std::array<std::array<LD, s>, s>;

    static constexpr arr  c = { ... };
    static constexpr arr  b = { ... };
    static constexpr arr  bstar = { ... };
    static constexpr LD   gamma = ...;
    static constexpr arr2 a = { ... };
    static constexpr arr2 g = { ... };
};
```

### Minimum coefficient requirements (Rosenbrock)

- `s`, `p`
- `a`, `g`
- `b`, `bstar`
- `c`
- scalar `gamma`

### Practical checklist (Rosenbrock)

1. Add method struct in `Rosenbrock/METHOD.hpp`.
2. Switch method macro in `Rosenbrock/Examples/makefile`, or use explicit solver alias in an example.
3. Build:
   ```bash
   make -C Rosenbrock/Examples
   ```
4. Run quick checks:
   ```bash
   ./Rosenbrock/Examples/Example.run | head
   ./Rosenbrock/Examples/Robertson.run | head
   ```

---

## 4) Common pitfalls when adding methods

1. **Wrong array sizes**  
   Ensure every coefficient array shape matches `s`.

2. **Inconsistent `c` values**  
   For many methods, `c[i]` should match row sums (where method definition expects that).

3. **Bad embedded pair (`b`/`bstar`)**  
   This directly affects adaptive error estimate quality.

4. **Numerical instability in stiff problems**  
   For Rosenbrock, verify `gamma`, `a`, `g` consistency and test against stiff benchmarks.

5. **Step controller too aggressive**  
   If new method is unstable, start with `simple` controller and conservative tolerances.

---

## 5) Hacking step control behavior

Step controller code lives here:

- RKF:
  - `RKF/RKF_step_control_simple.hpp`
  - `RKF/RKF_step_control_PI.hpp`
- Rosenbrock:
  - `Rosenbrock/Ros_step_control_simple.hpp`
  - `Rosenbrock/Ros_step_control_PI.hpp`

Selection is compile-time through:

```cpp
rkf::step_controllers::simple or ::PI
rosenbrock::step_controllers::simple or ::PI
```

If you add a new controller:

1. extend `step_controllers` enum,
2. add controller implementation function(s),
3. update the internal `step_control()` dispatcher (`if constexpr` chain),
4. test both existing and new paths.

---

## 6) Jacobian path (Rosenbrock-specific)

`rosenbrock::Solver` supports:

- default finite-difference Jacobian (constructor with `Jacobian_h`),
- user-provided Jacobian callback.

For performance-sensitive stiff systems, supplying analytic Jacobian callback usually gives better speed and robustness.

---

## 7) Suggested validation protocol for contributors

After solver/method changes:

```bash
make -C RKF/Examples
make -C Rosenbrock/Examples
make -C Example
```

And run:

```bash
./RKF/Examples/Example.run | head
./Rosenbrock/Examples/Example.run | head
./Rosenbrock/Examples/Robertson.run | head
./Example/TwoState.run 1 | head
```

When adding a new method, compare against at least one known solution/benchmark and inspect step count sensitivity to tolerances.

---

## 8) Documentation maintenance note

If you change solver signatures or add new public APIs:

1. update `docs/USAGE_AND_API.md`,
2. update this developer guide,
3. update examples if API usage changed.


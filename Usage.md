# NaBBODES

This document is a practical user guide.

---

## 1) How to use the solvers

NaBBODES currently provides two solver namespaces:

- `rkf` (explicit embedded Runge-Kutta)
- `rosenbrock` (Rosenbrock methods for stiff systems)

### 1.1 Include headers

```cpp
#include "RKF/RKF.hpp"

#include "Rosenbrock/Rosenbrock.hpp"
```

### 1.2 Define problem size and state type

```cpp
using LD = long double;
```

### 1.3 Provide RHS callback

Both solvers use the same RHS signature:

```cpp
void rhs(std::vector<LD>& lhs, const std::vector<LD>& y, const LD& t);
```

- `lhs` is output (`dy/dt`).
- `y` and `t` are input.

### 1.4 Choose a method and step controller

Step controllers:

- `step_controllers::simple`
- `step_controllers::PI`

Methods are defined in `RKF/METHOD.hpp` and `Rosenbrock/METHOD.hpp`.

### 1.5 Construct solver and run

```cpp
Solver system(rhs, y0, tmax, options...);
system.solve();
```

### 1.6 Read results

Use accessors:

- `get_t(...)`
- `get_solution(eq, step)`
- `get_error(eq, step)`
- `get_current_step()`

---

## 2) Public API reference

Below are all public member functions with signatures and descriptions.

---

## 2.1 `rkf::Solver`

Class template:

```cpp
template<class LD,, class RK_method = DormandPrince<LD>, 
          step_controllers step_controller = step_controllers::PI>
class Solver
```

### Public type aliases

```cpp
using system_type =
  std::function<void(std::vector<LD>& lhs,
                     const std::vector<LD>& y,
                     const LD& t)>;
```

### Constructors / destructor

```cpp
Solver(const system_type& dydt,
       const std::vector<LD>& init_cond,
       const LD& tmax,
       const parameters<LD>& opt = default_parameters<LD>);
```

Creates and initializes solver state with given RHS, initial condition, end time, and options. Note that the number of equations is determined from the size of init_cond.

```cpp
~Solver() = default;
```

### Integration control

```cpp
void next_step();
```
Advance by one accepted adaptive step.

```cpp
void solve();
```
Integrate from initial condition to `tmax`.

### Result accessors

```cpp
const std::vector<LD>& get_t() const;
```
Return full time grid.

```cpp
auto get_t(const unsigned int& step) const;
```
Return time at one step (bounds-checked via `.at()` in implementation).

```cpp
const std::vector<LD>& get_solution(const unsigned int& eq) const;
```
Return full trajectory of equation/component `eq`.

```cpp
auto get_solution(const unsigned int& eq, const unsigned int& step) const;
```
Return one solution value.

```cpp
const std::vector<LD>& get_error(const unsigned int& eq) const;
```
Return full embedded-error history for component `eq`.

```cpp
auto get_error(const unsigned int& eq, const unsigned int& step) const;
```
Return one embedded-error value.

### Parameter/state utilities

```cpp
void set_parameters(const parameters<LD>& opt = default_parameters<LD>);
```
Update runtime parameters (supports optional-field style updates).

```cpp
void reset(const std::vector<LD>& init_cond,
           const LD& tmax,
           const parameters<LD>& opt = default_parameters<LD>);
```
Reinitialize solver data and optionally update parameters.

```cpp
auto get_current_step() const;
```
Return current number of stored steps (size of time vector).

```cpp
auto get_current_step_size() const;
```
Return last accepted step size.

```cpp
const parameters<LD>& get_parameters() const;
```
Return active parameter set.

---

## 2.2 `rosenbrock::Solver`

Class template:

```cpp
template<class LD, class RK_method = RODAS5<LD>, 
          step_controllers step_controller = step_controllers::PI>
class Solver
```

### Public type aliases

```cpp
using system_type =
  std::function<void(std::vector<LD>& lhs,
                     const std::vector<LD>& y,
                     const LD& t)>;

using Jacobian_type =
  std::function<void(std::vector<std::vector<LD>>& J,
                     std::vector<LD>& dfdt,
                     const std::vector<LD>& y,
                     const LD& t)>;
```

### Constructors / destructor

Default finite-difference Jacobian:

```cpp
Solver(const system_type& dydt,
       const std::vector<LD>& init_cond,
       LD tmax,
       const parameters<LD>& opt = default_parameters<LD>,
       const LD& Jacobian_h = 1e-8);
```

Custom Jacobian callback:

```cpp
Solver(const system_type& dydt,
       const std::vector<LD>& init_cond,
       LD tmax,
       Jacobian_type Jac,
       const parameters<LD>& opt = default_parameters<LD>);
```

```cpp
~Solver() = default;
```

### Integration control

```cpp
void next_step();
void solve();
```

### Result accessors

```cpp
const std::vector<LD>& get_t() const;
auto get_t(const unsigned int& step) const;
const std::vector<LD>& get_solution(const unsigned int& eq) const;
auto get_solution(const unsigned int& eq, const unsigned int& step) const;
const std::vector<LD>& get_error(const unsigned int& eq) const;
auto get_error(const unsigned int& eq, const unsigned int& step) const;
```

### Parameter/state utilities

```cpp
void set_parameters(const parameters<LD>& opt = default_parameters<LD>);
void reset(const vector<LD>& init_cond,
           LD tmax,
           const parameters<LD>& opt = default_parameters<LD>);
auto get_current_step() const;
auto get_current_step_size() const;
const parameters<LD>& get_parameters() const;
```

---

## 3) Options (`parameters<LD>`)

Both solver namespaces define the same option fields (all `std::optional`):

- `initial_step_size`
- `minimum_step_size`
- `maximum_step_size`
- `maximum_No_steps`
- `absolute_tolerance`
- `relative_tolerance`
- `beta`
- `fac_max`
- `fac_min`

You can set only what you need:

```cpp
rkf::parameters<long double> opt{
    .absolute_tolerance = 1e-10L,
    .relative_tolerance = 1e-10L
};
```

---

## 4) Examples

## 4.1 Minimal RKF example

```cpp
#include <vector>
#include "RKF/RKF.hpp"

using LD = long double;

struct Eq {
    void operator()(std::vector<LD>& lhs, const std::vector<LD>& y, const LD& t) {
        lhs[0] = -y[0] + t;
    }
};

int main() {
    Eq f;
    Array y0{1.0};
    using S = rkf::Solver<LD, DormandPrince<LD>, rkf::step_controllers::PI>;
    S solver(f, y0, 1.0);
    solver.solve();
}
```

## 4.2 Minimal Rosenbrock example (default Jacobian)

```cpp
#include <vector>
#include "Rosenbrock/Rosenbrock.hpp"

using LD = long double;

struct Eq {
    void operator()(std::vector<LD>& lhs, const std::vector<LD>& y, const LD& t) {
        lhs[0] = -10.0L * y[0] + t;
    }
};

int main() {
    Eq f;
    std::vector y0(1,1.0);
    unsigned int N_eqs=y0.size();
    using S = rosenbrock::Solver<LD, ROS34PW2<LD>, rosenbrock::step_controllers::PI>;
    S solver(f, y0, 1.0); // default finite-difference Jacobian
    solver.solve();
}
```

## 4.3 Rosenbrock with custom Jacobian callback

```cpp
using S = rosenbrock::Solver<LD, ROS34PW2<LD>, rosenbrock::step_controllers::PI>;

auto jac = [](std::vector<vector<LD>>& J, std::vector<LD>& dfdt, const vector<LD>& y, const LD& t) {
    J.assign(N_eqs,std::vector<LD>(N_eqs,-10));
    dfdt[0].assign(N_eqs,std::vector<LD>(N_eqs,-1)); // d/dt f(y,t)
};

S solver(f, y0, 1.0, jac);
solver.solve();
```

## 4.4 Build and run shipped examples

```bash
make -C RKF/Examples
make -C Rosenbrock/Examples
make -C Example

./RKF/Examples/Example.run
./Rosenbrock/Examples/Example.run
./Example/TwoState.run 1
```
---

# Adding new methods


## Add a new RKF method 

Different methods are available in`RKF/METHOD.hpp`. New methods can be defined a struct:

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

The variables should be given as `constexpr` because they are not meant to change.

The different variables correspond to how RKF takes the iterative steps:

Assuming that we have a system of differential equations

$$
\dfrac{d\vec{y}}{dt}=\vec{f}(\vec{y},t) \;,
$$

with  given $\vec{y}\left(0\right)$ (assume that integration starts at $t=0$). 

RKF follows the iteration

$$
\vec{y}_{n+1}=\vec{y}_{n}+ h\sum_{i=1}^{s} b_i \vec{k}_i
$$
$$
\vec{y}^{\star}_{n+1}=\vec{y}_{n}+ h\sum_{i=1}^{s} b_i^{\star} \vec{k}_i \;,
$$
with
$$
\vec{k}_{i}=\vec f\Bigg(\vec{y}_{n}+h \Big(\sum_{j=1}^{i-1}a_{ij}\vec{k}_{j} \Big), t_{n}+h c_{i}   \Bigg)\;.
$$




## Add a new Rosenbrock method 

Different methods are available in`Rosenbrock/METHOD.hpp`. New methods can be defined a struct:

```cpp
template<class LD>
struct MyRosMethod {
    static constexpr unsigned int s = ...;   // stage
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

The variables should be given as `constexpr` because they are not meant to change.

The different variables correspond to how Rosenbrock takes the iterative steps:


Assuming that we have a system of differential equations

$$
\dfrac{d\vec{y}}{dt}=\vec{f}(\vec{y},t) \;,
$$

with given $\vec{y}\left(0\right)$ (assume that integration starts at $t=0$). 


After a bit ao algebra, the Rosenbrock method solves differential equations using 

$$
\vec{y}_{n+1}=\vec{y}_{n}+ \sum_{i=1}^{s} b_i \vec{k}_i
$$
$$
\vec{y}^{\star}_{n+1}=\vec{y}_{n}+ \sum_{i=1}^{s} b_i^{\star} \vec{k}_i 
$$

$$
\left(\hat I - \gamma h \hat J\right)\cdot \vec{k}_{i}=
h \vec{f}\Big(\vec{y}_n+\sum_{j=1}^{i-1}a_{ij}\vec{k}_{j},t_n + c_{i}h \Big)+
h^2 \left(\gamma + \sum_{j=1}^{i-1}\gamma_{ij}\right)\dfrac{\partial \vec{f} }{\partial t}+
h \hat J \cdot \sum_{j=1}^{i-1}\gamma_{ij} \vec{k}_{j}\; ,
$$

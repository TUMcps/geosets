| Operation                              | HPolytope | VPolytope | Polygon | Interval | Zonotope |
| -------------------------------------- | --------- | --------- | ------- | -------- | -------- |
| Empty check                            | ✔         | ✔         | ✔       | ✔        | ✔        |
| Bounded check                          | ✔         | ✔         | ✔       | ✔        | ✔        |
| Center computation                     | ✔         | ✔         | ✔       | ✔        | ✔        |
| Volume                                 | ✔**       | ✔         | ✔       | ✔        | ✔        |
| Containment check                      | ✔         | ✔         | ✔       | ✔        | ✔        |
| Boundary point                         | ✔         | ✔         | :x:     | ✔        | ✔        |
| Boundary points from line within set   | ✔         | ✔         | ✔       | ✔        | ✔        |
| Conversion to vertices (plotting)      | ✔         | ✔         | ✔       | ✔        | ✔        |
| Translation w. vector                  | ✔         | ✔         | ✔       | ✔        | ✔        |
| Linear transform                       | ✔         | ✔         | ✔       | ✔        | ✔        |
| Support function                       | ✔         | ✔         | ✔       | ✔        | ✔        |
| Intersect predicate (w. same set rep.) | ✔         | ✔**       | ✔       | ✔        | :x:      |
| Intersection (w. same set rep.)        | ✔         | ✔**       | ✔       | ✔        | :x:      |
| Minkowski sum (w. same set rep.)       | ✔         | ✔         | ✔       | ✔        | ✔        |

\* only for specific dimensions
** uses conversion to another representation

### Conversion

| Conversion     | HP  | VP  | PG  | I   | Z   |
| -------------- | --- | --- | --- | --- | --- |
| HPolytope (HP) |     | ✔   | :x: | :x: | ✔   |
| VPolytope (VP) | ✔   |     | ✔   | :x: | :x: |
| Polygon (PG)   | :x: | ✔   |     | :x: | :x: |
| Interval (I)   | :x: | :x: | :x: |     | :x: |
| Zonotope (Z)   | :x: | :x: | :x: | :x: |     |

### Sampling

| Operation                          | HPolytope | VPolytope | Polygon | Interval | Zonotope |
| ---------------------------------- | --------- | --------- | ------- | -------- | -------- |
| Gaussian RDHR sampling             | ✔         | :x:       | ✔       | ✔        | :x:      |
| Gaussian rejection sampling        | ✔         | ✔         | ✔       | ✔        | ✔        |
| Gaussian rejection + RDHR fallback | ✔         | :x:       | ✔       | ✔        | :x:      |

### Representation-specific functions:

Polygon:

- set difference

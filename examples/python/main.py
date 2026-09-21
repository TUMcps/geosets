import numpy as np

import geosets_py as gsp


def operate_on_set(set_obj: gsp.GeometricSet):
    gsp.print_name(set_obj)

    print(set_obj.dim())
    print(set_obj.translate(np.zeros(2)))
    print(set_obj.translate_(np.zeros(2)))


def main():
    hpoly = gsp.HPolytope.from_random(2, 4)
    vpoly = gsp.VPolytope.from_random(2, 10)
    polygon = gsp.Polygon.from_random(2, 10)
    interval = gsp.Interval.from_random(2)
    _ = interval.translate(np.zeros(2))

    vpoly.plot(fill_opacity=0.1)
    interval.plot()
    polygon.plot(show=True, color="green")

    print(type(hpoly))
    print(isinstance(hpoly, gsp.HPolytope))
    print(isinstance(hpoly, gsp.GeometricSet))
    print(type(hpoly) is gsp.GeometricSet)

    print("Function calling")
    operate_on_set(hpoly)
    operate_on_set(vpoly)


if __name__ == "__main__":
    main()

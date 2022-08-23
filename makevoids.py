import matplotlib.pyplot as plt
import numpy as np
from collections import namedtuple

# Some immutable structs
Rectangle = namedtuple("Domain", ["xmin", "ymin", "xmax", "ymax"])
Void_stats = namedtuple("Void_stats", ["n_voids", "rad_mean", "rad_std", "min_spacing"])
Void = namedtuple("Void", ["x", "y", "radius"])

# Use a quasi-random generator for reproducability
rng = np.random.default_rng(12341)


def gen_void_from_stats(void_stats, bounds):
    domain_width = bounds.xmax - bounds.xmin
    domain_heigth = bounds.ymax - bounds.ymin

    diameter = np.abs(rng.normal(loc=void_stats.rad_mean, scale=void_stats.rad_std))
    x = rng.random() * domain_width + bounds.xmin
    y = rng.random() * domain_heigth + bounds.ymin

    return Void(x, y, diameter)


def calc_void_spacing(void1, void2):
    return (
        np.sqrt((void1.x - void2.x) ** 2.0 + (void1.y - void2.y) ** 2.0)
        - void1.radius
        - void2.radius
    )


def void_is_to_close(void, voids, min_spacing):
    if len(voids) == 0:
        # No possible collisions
        return False
    else:
        for existing_void in voids:
            if calc_void_spacing(void, existing_void) < min_spacing:
                return True
        return False


def void_is_within_bounds(void, domain_shape):
    void_x_min = void.x - void.radius
    void_x_max = void.x + void.radius
    void_y_min = void.y - void.radius
    void_y_max = void.y + void.radius

    if (
        void_x_min < domain_shape.xmin
        or void_x_max > domain_shape.xmax
        or void_y_min < domain_shape.ymin
        or void_y_max > domain_shape.ymax
    ):
        return False
    else:
        return True


def make_voided_domain(void_stats, bounds, max_attempts=50):
    voids = []
    for _ in range(void_stats.n_voids):
        for _ in range(max_attempts):
            void = gen_void_from_stats(void_stats, bounds)
            if not void_is_to_close(
                void, voids, min_spacing=void_stats.min_spacing
            ) and void_is_within_bounds(void, bounds):
                voids.append(void)
                break
    return voids


def plot_voids(voids):
    figure, ax = plt.subplots()
    for void in voids:
        circle = plt.Circle((void.x, void.y), void.radius)
        ax.add_artist(circle)

    ax.set_title("Voided domain")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_aspect("equal")
    plt.xlim(left=bounds.xmin, right=bounds.xmax)
    plt.ylim(top=bounds.ymin, bottom=bounds.ymax)
    plt.show()


# Rectangular domain where the voids are to be inserted
bounds = Rectangle(xmin=250, ymin=30, xmax=750, ymax=570)
# bounds = Rectangle(xmin=20, ymin=30, xmax=970, ymax=580)
# Statistics of the voids to be inserted
# void_stats = Void_stats(n_voids=700, rad_mean=2, rad_std=0.0, min_spacing=3)
void_stats = Void_stats(n_voids=10, rad_mean=20, rad_std=0.0, min_spacing=3)

voids = make_voided_domain(void_stats, bounds)
print(
    "Void volume fraction is: ",
    np.sum([np.pi * void.radius ** 2 for void in voids])
    / ((bounds.xmax - bounds.xmin) * (bounds.ymax - bounds.ymin)),
)

plot_voids(voids)

np.savetxt("voids.csv", voids, delimiter=",", header="x,y,diameter")


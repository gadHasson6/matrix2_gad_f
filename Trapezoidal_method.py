
from colors import bcolors



def trapezoidal_rule(f, a, b, n):
    if a > b:
        a, b = b, a
    h = (b - a) / n
    T = f(a) + f(b)
    integral = ((b-a) * T) / 2  # Initialize with endpoints

    for i in range(1, n):
        x_i = a + i * h
        integral += f(x_i)

    integral *= h

    return abs(integral)


if __name__ == '__main__':
    f = lambda x: x**2
    result = trapezoidal_rule(f, 0, 1, 10)
    print(bcolors.OKBLUE, "Approximate integral:", result, bcolors.ENDC)

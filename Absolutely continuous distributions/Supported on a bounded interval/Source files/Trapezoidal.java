public class Trapezoidal extends MathLibrary{
    double a;
    double b;
    double c;
    double d;

    Trapezoidal(double a, double b, double c, double d) {
        if (a >= d) throw new IllegalArgumentException("a < d");
        if (b < a || b >= c) throw new IllegalArgumentException("a <= b < c");
        if (c <= b || c > d) throw new IllegalArgumentException("b < c <= d");
        if (d > c) throw new IllegalArgumentException("c <= d");
        this.a = a;
        this.b = b;
        this.c = c;
        this.d = d;
    }

    double PDF(double x) {
        if (x < a || x > d) throw new IllegalArgumentException("a <= x <= b");
        if (a <= x && x < b) {
            return (2 / (d + c - a - b)) * ((x - a) / (b - a));
        } else if (b <= x && x < c) {
            return 2 / (d + c - a - b);
        } else {
            return (2 / (d + c - a - b)) / ((d - x) / (d - c));
        }
    }

    double CDF(double x) {
        if (x < a || x > d) throw new IllegalArgumentException("a <= x <= b");
        double part = 1 / (d + c - a - b);
        if (a <= x && x < b) {
            return part * (1 / (b - a)) * power(x - a, 2);
        } else if (b <= x && x < c) {
            return part * (2 * x - a - b);
        } else {
            return 1 - (part * (1 / (d - c)) * power(d - x, 2));
        }
    }

    double Mean() {
        double part1 = 1 / (3 * (d + c - b - a));
        double part2 = ((power(d, 3) - power(c, 3)) / (d - c)) - ((power(b, 3) - power(a, 3)) / (b - a));
        return part1 * part2;
    }

    double Variance() {
        double part1 = 1 / (6 * (d + c - b - a));
        double part2 = ((power(d, 4) - power(c, 4)) / (d - c)) - ((power(b, 4) - power(a, 4)) / (b - a));
        return part1 * part2;
    }

    double Entropy() {
        double part1 = (d - c + b - a) / (2 * (d + c - b - a));
        double part2 = ln((d + c - b - a) / 2);
        return part1 + part2;
    }

    double MGF(double t) {
        double part1 = (2 / (d + c - b - a)) * (1 / power(t, 2));
        double part2 = ((euler(d * t) - euler(c * t)) / (d - c)) - ((euler(b * t) - euler(a * t)) / (b - a));
        return part1 * part2;
    }
}

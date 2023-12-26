public class Triangular extends MathLibrary{
    double a;
    double b;
    double c;

    Triangular(double a, double b, double c) {
        if (b <= a) throw new IllegalArgumentException(" a < b");
        if (c < a || c > b) throw new IllegalArgumentException("a <= c <= b");
        this.a = a;
        this.b = b;
        this.c = c;
    }

    double PDF(double x) {
        if (x < a) {
            return 0;
        } else if (x >= a && x < c) {
            return (2 * (x - a)) / ((b - a) * (c - a));
        } else if (x == c) {
            return 2 / (b - a);
        } else if (x > c && x <= b) {
            return (2 * (b - x)) / ((b - a) * (b - c));
        } else {
            return 0;
        }
    }

    double CDF(double x) {
        if (x <= a) {
            return 0;
        } else if (x > a && x <= c) {
            return power(x - a, 2) / ((b - a) * (c - a));
        } else if (x > c && x < b) {
            return 1 - (power(b - a, 2) / ((b - a) * (b - c)));
        } else {
            return 1;
        }
    }

    double Mean() {
        return (a + b + c) / 3;
    }

    double Median() {
        if (c >= ((a + b) / 2)) {
            return a + root(((b - a) * (c - a)) / 2, 2);
        } else {
            return b - root(((b - a) * (b - c)) / 2, 2);
        }
    }

    double Mode() {
        return c;
    }

    double Variance() {
        return (power(a, 2) + power(b, 2) + power(c, 2) - (a * b) - (a * c) - (b * c));
    }

    double Skewness() {
        double part1 = root(2, 2) * (a + b - (2 * c)) * (2 * a - b - c) * (a - (2 * b) + c);
        double part2x = power(a, 2) + power(b, 2) + power(c, 2) - (a * b) - (a * c) - (b * c);
        double part2 = 5 * part2x * root(part2x, 2);
        return part1 / part2;
    }

    double EXkurtosis() {
        return -0.6;
    }

    double Entropy() {
        return 0.5 + ln((b - a) / 2);
    }

    double MGF(double t) {
        double part1 = ((b - c) * euler(a * t)) - ((b - a) * euler(c * t)) + ((c - a) * euler(b * t));
        double part2 = (b - a) * (c - a) * (b - c) * power(t, 2);
        return 2 * part1 / part2;
    }
}

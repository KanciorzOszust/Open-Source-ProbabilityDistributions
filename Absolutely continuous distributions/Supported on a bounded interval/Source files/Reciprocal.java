public class Reciprocal extends MathLibrary{
    double a;
    double b;

    Reciprocal(double a, double b) {
        if (a <= 0 || a >= b) throw new IllegalArgumentException("0 < a < b");
        this.a = a;
        this.b = b;
    }

    double PDF(double x) {
        if (x < a || x > b) throw new IllegalArgumentException("a <= x <= b");
        return 1 / (x * ln(b / a));
    }

    double CDF(double x) {
        if (x < a || x > b) throw new IllegalArgumentException("a <= x <= b");
        return ln(x / a) / ln(b / a);
    }

    double Mean() {
        return (b - a) / ln(b / a);
    }

    double Median() {
        return root(a * b, 2);
    }

    double Variance() {
        double part1 = (power(b, 2) - power(a, 2)) / (2 * ln(b / a));
        double part2 = power((b - a) / ln(b / a), 2);
        return part1 - part2;
    }
}

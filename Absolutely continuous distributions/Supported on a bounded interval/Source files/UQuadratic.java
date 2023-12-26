public class UQuadratic extends MathLibrary{
    double a;
    double b;

    UQuadratic(double a, double b) {
        if (b <= a) throw new IllegalArgumentException("b > a");
        this.a = a;
        this.b = b;
    }

    double PDF(double x) {
        if (x < a || x > b) throw new IllegalArgumentException("a <= x <= b");
        return a * power(x - b, 2);
    }

    double CDF(double x) {
        if (x < a || x > b) throw new IllegalArgumentException("a<= x <= b");
        return (a / 3) * (power(x - b, 3) + power(b - a, 3));
    }

    double Mean() {
        return (a + b) / 2;
    }

    double Median() {
        return Mean();
    }

    double[] Mode() {
        double[] arr = {a, b};
        return arr;
    }

    double Variance() {
        return 0.15 * power(b - a, 2);
    }

    int Skewness() {
        return 0;
    }

    double EXkurtosis() {
        return (3 / 112.0) * power(b - a, 4);
    }

    double MGF(double t) {
        double part1 = -3 * (euler(a * t) * (4 + (power(a, 2) + (2 * a) * (-2 + b) + power(b, 2)) * t));
        double part2 = euler(b * t) * (4 + (-4 * b + power(a + b, 2)) * t);
        return (part1 - part2) / (power(a - b, 3) * power(t, 2));
    }
}

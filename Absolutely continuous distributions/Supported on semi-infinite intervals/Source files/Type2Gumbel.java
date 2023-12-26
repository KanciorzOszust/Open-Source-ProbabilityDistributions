public class Type2Gumbel extends MathLibrary{
    double a;
    double b;

    Type2Gumbel(double a, double b) {
        this.a = a;
        this.b = b;
    }

    double PDF(double x) {
        if (x <= 0) throw new IllegalArgumentException("x > 0");
        return a * b * doublePower(x, -1 * a - 1) * euler(-1 * b * doublePower(x, -1 * a));
    }

    double CDF(double x) {
        if (x <= 0) throw new IllegalArgumentException("x > 0");
        return euler(-1 * b * doublePower(x, -1 * a));
    }

    double Mean() {
        return doublePower(b, 1 / a) * gamma(1 - (1 / a));
    }

    double Variance() {
        return doublePower(b, 2 / a) * (gamma(1 - (1 / a) - power(gamma(1 - (1 / a)), 2)));
    }
}

public class Kumaraswamy extends MathLibrary{
    double a;
    double b;

    Kumaraswamy(double a, double b) {
        if (a <= 0 || b <= 0) throw new IllegalArgumentException("a,b > 0");
        this.a = a;
        this.b = b;
    }

    double PDF(double x) {
        if (x <= 0 || x >= 1) throw new IllegalArgumentException("0 < x < 1");
        return a * b * doublePower(x, a - 1) * doublePower(1 - doublePower(x, a), b - 1);
    }

    double CDF(double x) {
        if (x <= 0 || x >= 1) throw new IllegalArgumentException("0 < x < 1");
        return 1 - doublePower(1 - doublePower(x, a), b - 1);
    }

    double Mean() {
        return (b * gamma(1 + (1 / a)) * gamma(b)) / gamma(1 + (1 / a) + b);
    }

    double Median() {
        return doublePower(1 - doublePower(0.5, (1 / b)), (1 / a));
    }

    double Mode() {
        if (a <= 1 || b <= 1) throw new IllegalArgumentException("a, b > 1");
        return doublePower((a - 1) / (a * b - 1), (1 / a));
    }

    double Entropy() {
        return (1 - (1/ b)) + ((1 - (1 / a)) * generalizedHarmonicNumber(1, b)) - ln(a * b);
    }
}

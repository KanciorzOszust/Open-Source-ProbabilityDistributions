public class GeneralizedGamma extends MathLibrary{
    double a;
    double d;
    double p;

    GeneralizedGamma(double a, double d, double p) {
        if (a <= 0 || d <= 0 || p <= 0) throw new IllegalArgumentException("a, d, p > 0");
        this.a = a;
        this.d = d;
        this.p = p;
    }

    double PDF(double x) {
        if (x <= 0) throw new IllegalArgumentException("x > 0");
        double part1 = p / doublePower(a, d) / gamma(d / p);
        double part2 = doublePower(x, d - 1) * euler(-1 * doublePower(x / a, p));
        return part1 * part2;
    }

    double CDF(double x) {
        if (x <= 0) throw new IllegalArgumentException("x > 0");
        return incompleteLowerGamma(d / p, doublePower(x / a, p)) / gamma(d / p);
    }

    double Mean() {
        return a * gamma((d + 1) / p) / (gamma(d / p));
    }

    double Mode() {
        if (d > 1) {
            return a * doublePower((d - 1) / p, 1 / p);
        } else {
            return 0;
        }
    }

    double Variance() {
        double part1 = gamma((d + 2) / p) / gamma(d / p);
        double part2 = power(gamma((d + 1) / p) / gamma(d / p), 2);
        return power(a, 2) * (part1 - part2);
    }

    double Entropy() {
        return ln((a * gamma(d / p)) / p) + (d / p) + (((1 / p) - (d / p)) * digamma(d / p));
    }
}

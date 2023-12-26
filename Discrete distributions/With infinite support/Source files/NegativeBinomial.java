public class NegativeBinomial extends MathLibrary{
    int r;
    double p;

    NegativeBinomial(int r, double p) {
        if (r <= 0) throw new IllegalArgumentException("r > 0");
        if (p < 0 || p > 1) throw new IllegalArgumentException("0 <= p <= 1");
        this.r = r;
        this.p = p;
    }

    double PMF(int k) {
        return binomial(k + r - 1, k) * power(1 - p, k) * power(p, r);
    }

    double CDF(int k) {
        return regularizedIncompleteBeta(r, k + 1, p);
    }

    double Mean() {
        return (r * (1 - p)) / p;
    }

    int Mode() {
        if (r > 1) {
            return (int) (((r - 1) * (1 - p)) / p);
        } else {
            return 0;
        }
    }

    double Variance() {
        return (r * (1 - p)) / power(p, 2);
    }

    double Skewness() {
        return (2 - p) / root((1 - p) * r, 2);
    }

    double EXkurtosis() {
        return 6.0 / r + (power(p, 2) / ((1 - p) * r));
    }

    double MGF(double t) {
        if (t < (-1 * log(1 - p, 10))) {
            double base = p / (1 - ((1 - p) * euler(t)));
            return power(base, r);
        } else {
            return 0;
        }
    }

    double PGF(double z) {
        if (absolute(z) < (1 / p)) {
            double base = p / (1 - ((1 - p) * z));
            return power(base, r);
        } else {
            return 0;
        }
    }

    double FisherInformation() {
        return r / (power(p, 2) * (1 - p));
    }
}

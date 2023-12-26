public class YuleSimon extends MathLibrary{
    double p;

    YuleSimon(double p) {
        if (p <= 0) throw new IllegalArgumentException("p > 0");
        this.p = p;
    }

    double PMF(int k) {
        return p * beta(k, p + 1);
    }

    double CDF(int k) {
        return 1 - (k * beta(k, p + 1));
    }

    double Mean() {
        if (p > 1) {
            return p / (p - 1);
        } else {
            return 0;
        }
    }

    double Mode() {
        return 1;
    }

    double Variance() {
        if (p > 2) {
            return power(p, 2) / (power(p - 1, 2) * (p - 2));
        } else {
            return 0;
        }
    }

    double Skewness() {
        if (p > 3) {
            return (power(p + 1, 2) * root(p - 2, 2)) / ((p - 3) * p);
        } else {
            return 0;
        }
    }

    double EXkurtosis() {
        if (p > 4) {
            return p + 3 + (((11 * power(p, 3)) - (49 * p) - 22) / ((p - 4) * (p - 3) * p));
        } else {
            return 0;
        }
    }
}

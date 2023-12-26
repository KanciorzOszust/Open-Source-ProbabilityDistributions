public class Stable extends MathLibrary{
    double alpha;
    double beta;
    double c;
    double mu;

    Stable(double alpha, double beta, double c, double mu) {
        if (alpha <= 0 || alpha > 2) throw new IllegalArgumentException("0 < alpha <= 2");
        if (beta < -1 || beta > 1) throw new IllegalArgumentException("-1 <= beta <= 1");
        if (c <= 0) throw new IllegalArgumentException("c > 0");
        this.alpha = alpha;
        this.beta = beta;
        this.c = c;
        this.mu = mu;
    }

    double Mean() {
        if (alpha > 1) return mu;
        else return Double.NaN;
    }

    double Median() {
        if (beta == 0) return mu;
        else return Double.NaN;
    }

    double Mode() {
        return Median();
    }

    double Variance() {
        if (alpha == 2) return 2 * power(c, 2);
        else return Double.NaN;
    }

    double Skewness() {
        if (alpha == 2) return 0;
        else return Double.NaN;
    }

    double EXkurtosis() {
        return Skewness();
    }

    double MGF(double t) {
        if (alpha == 2) {
            return euler(t * mu + (power(c, 2) * power(t, 2)));
        } else {
            return Double.NaN;
        }
    }
}

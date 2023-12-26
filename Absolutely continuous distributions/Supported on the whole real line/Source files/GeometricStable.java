public class GeometricStable extends MathLibrary{
    double alpha;
    double beta;
    double lambda;
    double mu;

    GeometricStable(double alpha, double beta, double lambda, double mu) {
        if (alpha <= 0 || alpha > 2) throw new IllegalArgumentException("0 < alpha <= 2");
        if (beta < -1 || beta > 1) throw new IllegalArgumentException("-1 <= beta <= 1");
        if (lambda <= 0) throw new IllegalArgumentException("lambda > 0");
        this.alpha = alpha;
        this.beta = beta;
        this.lambda = lambda;
        this.mu = mu;
    }

    double Median() {
        if (beta == 0) return mu;
        else return Double.NaN;
    }

    double Mode() {
        return Median();
    }

    double Variance() {
        if (alpha == 2) return 2 * power(lambda, 2);
        else return Double.NaN;
    }

    double Skewness() {
        if (alpha == 2) return 0;
        else return Double.NaN;
    }

    double EXkurtosis() {
        if (alpha == 2) return 3;
        else return Double.NaN;
    }
}

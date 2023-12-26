public class NormalExponentialGamma extends MathLibrary{
    double mu;
    double k;
    double theta;

    NormalExponentialGamma(double mu, double k, double theta) {
        if (k <= 0 || theta <= 0) throw new IllegalArgumentException("k, theta > 0");
        this.mu = mu;
        this.k = k;
        this.theta = theta;
    }

    double Mean() {
        return mu;
    }

    double Median() {
        return mu;
    }

    double Mode() {
        return mu;
    }

    double Variance() {
        return power(theta, 2) / (k - 1);
    }

    int Skewness() {
        return 0;
    }
}

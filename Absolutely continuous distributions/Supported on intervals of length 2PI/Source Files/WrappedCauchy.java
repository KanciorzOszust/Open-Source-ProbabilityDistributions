public class WrappedCauchy extends MathLibrary{
    double mu;
    double gamma;

    WrappedCauchy(double mu, double gamma) {
        if (gamma <= 0) throw new IllegalArgumentException("gamma > 0");
        this.mu = mu;
        this.gamma = gamma;
    }

    double PDF(double x) {
        if (x < -1 * PI || x >= PI) throw new IllegalArgumentException("-PI <= x < PI");
        return (1 / (2 * PI)) * (hyperbolicSine(gamma) / (hyperbolicCosine(gamma) - cosine(x - mu)));
    }

    double Mean() {
        return mu;
    }

    double Variance() {
        return 1 - euler(-1 * gamma);
    }

    double Entropy() {
        return ln(2 * PI * (1 - euler(-2 * gamma)));
    }
}
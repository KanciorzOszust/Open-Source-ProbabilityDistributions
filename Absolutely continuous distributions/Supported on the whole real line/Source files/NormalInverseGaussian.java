public class NormalInverseGaussian extends MathLibrary{
    double mu;
    double alpha;
    double beta;
    double delta;
    double gamma;

    NormalInverseGaussian(double mu, double alpha, double beta, double delta) {
        this.mu = mu;
        this.alpha = alpha;
        this.beta = beta;
        this.delta = delta;
        this.gamma = root(power(alpha, 2) - power(beta, 2), 2);
    }

    double PDF(double x) {
        double part1 = alpha * delta * modifiedBessel2(1, alpha * root(power(delta, 2) + power(x - mu, 2), 2));
        double part2 = PI * root(power(delta, 2) + power(x - mu, 2), 2);
        double part3 = euler(delta * gamma + (beta * (x - mu)));
        return part1 / part2 * part3;
    }
    
    double Mean() {
        return mu + (delta * beta / gamma);
    }

    double Variance() {
        return delta * power(alpha, 2) / power(gamma, 3);
    }

    double Skewness() {
        return (3 * beta) / (alpha * root(delta * gamma, 2));
    }

    double EXkurtosis() {
        return (3 * (1 + ((4 * power(beta, 2)) / power(alpha, 2)))) / (delta * gamma);
    }

    double MGF(double t) {
        return euler(mu * t + (delta * (gamma - root(power(alpha, 2) - power(beta + t, 2), 2))));
    }
}

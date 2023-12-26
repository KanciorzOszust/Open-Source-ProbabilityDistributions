public class Hyperbolic extends MathLibrary{
    double mu;
    double alpha;
    double beta;
    double delta;
    double gamma;

    Hyperbolic(double mu, double alpha, double beta, double delta) {
        this.mu = mu;
        this.alpha = alpha;
        this.beta = beta;
        this.delta = delta;
        this.gamma = root(power(alpha, 2) - power(beta, 2), 2);
    }

    double PDF(double x) {
        double part1 = gamma / (2 * alpha * delta * modifiedBessel2(1, delta * gamma));
        double part2 = -1 * alpha * root(power(delta, 2) + power(x - mu, 2), 2) + (beta * (x - mu));
        return part1 * euler(part2);
    }

    double Mean() {
        return mu + ((delta * beta * modifiedBessel2(2, delta * gamma)) / (gamma * modifiedBessel2(1, delta * gamma)));
    }

    double Mode() {
        return mu + ((delta * beta) / gamma);
    }

    double Variance() {
        double part1 = (delta * modifiedBessel2(2, delta * gamma)) / (gamma * modifiedBessel2(1, delta * gamma));
        double part2 = (power(beta, 2) * power(delta, 2)) / power(gamma, 2);
        double part3 = modifiedBessel2(3, delta * gamma) / modifiedBessel2(1, delta * gamma);
        double part4 = power(modifiedBessel2(2, delta * gamma), 2) / power(modifiedBessel2(1, delta * gamma), 2);
        return part1 + (part2 * (part3 - part4));
    }

    double MGF(double t) {
        double part1 = euler(mu * t) * gamma * modifiedBessel2(1, delta * root(power(alpha, 2) - power(beta + t, 2), 2));
        double part2 = root(power(alpha, 2) - power(beta + t, 2), 2) * modifiedBessel2(1, delta * gamma);
        return part1 / part2;
    }
}

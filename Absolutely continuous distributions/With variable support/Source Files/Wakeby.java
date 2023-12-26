public class Wakeby extends MathLibrary{
    double alpha;
    double beta;
    double gamma;
    double delta;
    double xi;

    Wakeby(double alpha, double beta, double gamma, double delta, double xi) {
        this.alpha = alpha;
        this.beta = beta;
        this.gamma = gamma;
        this.delta = delta;
        this.xi = xi;
    }

    double Quantile(double p) {
        if (delta >= 0 && gamma > 0) {
            if (p < xi) throw new IllegalArgumentException("p >= xi");
        } else {
            if (p < xi || p > (xi + (alpha / beta) - (gamma / delta))) throw new IllegalArgumentException("xi <= p <= xi + (alpha / beta) - (gamma / delta)");
        }
        double q = 1 - p;
        double part1 = xi + ((alpha / beta) * (1 - doublePower(q, beta)));
        double part2 = (gamma / delta) * (1 - doublePower(q, -1 * delta));
        return part1 - part2;
    }
}

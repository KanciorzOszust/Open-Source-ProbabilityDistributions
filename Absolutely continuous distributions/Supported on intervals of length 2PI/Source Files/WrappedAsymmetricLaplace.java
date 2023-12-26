public class WrappedAsymmetricLaplace extends MathLibrary{
    double m;
    double gamma;
    double k;
    
    WrappedAsymmetricLaplace(double m, double gamma, double k) {
        if (m < 0 || m >= 2 * PI) throw new IllegalArgumentException("0 <= m < 2PI");
        if (gamma <= 0) throw new IllegalArgumentException("gamma > 0");
        if (k <= 0) throw new IllegalArgumentException("k > 0");
        this.m = m;
        this.gamma = gamma;
        this.k = k;
    }

    double PDF(double x) {
        if (m < 0 || m >= 2 * PI) throw new IllegalArgumentException("0 <= x < 2PI");
        double standard = (k * gamma) / (power(k, 2) + 1);
        double part1;
        double part2;
        if (x >= m) {
            part1 = euler(-1 * (x - m) * gamma * k) / (1 - euler(-2 * PI * gamma * k));
            part2 = euler(-1 * (x - m) * (gamma / k)) / (1 - euler(2 * PI * (gamma / k)));
        } else {
            part1 = euler(-1 * (x - m) * gamma * k) / (euler(2 * PI * gamma * k) - 1);
            part2 = euler((x - m) * gamma / k) / (euler(-2 * PI * gamma / k) - 1);
        }
        return standard * (part1 - part2);
    }

    double Mean() {
        return m;
    }

    double Variance() {
        double part = (power(k, -2) + power(gamma, 2)) * (power(k, 2) + power(gamma, 2));
        return 1 - (power(gamma, 2) / power(part, 2));
    }
}

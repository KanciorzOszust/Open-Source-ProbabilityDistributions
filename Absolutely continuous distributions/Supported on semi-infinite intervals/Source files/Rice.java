public class Rice extends MathLibrary{
    double nu;
    double sigma;
    double sigmaSquare;

    Rice(double nu, double sigma) {
        if (nu < 0 || sigma < 0) throw new IllegalArgumentException("nu, sigma >= 0"); 
        this.nu = nu;
        this.sigma = sigma;
        this.sigmaSquare = power(sigma, 2);
    }

    double PDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        return (x / sigmaSquare) * euler((-1 * (power(x, 2) + power(nu, 2))) / (2 * sigmaSquare)) * modifiedBessel(0, x * nu / sigmaSquare);
    }

    double CDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        return 1 - marcumQ(1, nu / sigma, x / sigma);
    }

    double Mean() {
        return (sigma * root(PI / 2, 2)) * laguerre(0.5, 1, -1 * power(nu, 2) / (2 * sigmaSquare));
    }

    double Variance() {
        double part = (PI * sigmaSquare / 2) * laguerre(0.5, 2, -1 * power(nu, 2) / (2 * sigmaSquare));
        return 2 * sigmaSquare + power(nu, 2) - part;
    }

}

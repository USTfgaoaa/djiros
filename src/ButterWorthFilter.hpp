#include <cstdlib>
#include <cassert>
#include <vector>
#include <deque>
#include <cmath>
#include <limits>
#include <stdexcept>

class ButterWorthFilter {
  public:
    // Filter parameters
    const int order = 2;
    const double sample_rate = 100.0;

    // H(z) = B(z) / A(z) = [b(0) + b(1)z^-1 + ... + b(n)z^-n] / [a(0) + a(1)z^-1 + ... + a(n)z^-n]
    const std::vector<double> a = {1.000000000000000, -1.561018075800718, 0.641351538057563};
    const std::vector<double> b = {0.020083365564211, 0.040166731128423, 0.020083365564211};

    // Group delay in the passband. Delay in the number of samples at sample_rate.
    const int avg_delay = 5; 

  public:
    ButterWorthFilter() : sample_rate_verified(false){};

    bool update(double x, double& y) {
        if (!sample_rate_verified) {
            throw std::runtime_error("Sample rate is not verified!");
        }

        x_buf.push_front(x);

        if (x_buf.size() <= static_cast<size_t>(order)) {  // not enough data
            y_buf.push_front(0.0);
            y = std::numeric_limits<double>::quiet_NaN();
            return false;
        }

        // Get enough data
        y_buf.push_front(b[0] * x_buf[0]);
        for (int i = 1; i <= order; ++i) {
            y_buf[0] = y_buf[0] + b[i] * x_buf[i] - a[i] * y_buf[i];
        }
        y = y_buf[0];
        // std::cout << y_buf[0] << " " << y_buf[1] << " " << y_buf[2] << std::endl;

        x_buf.pop_back();
        y_buf.pop_back();
        assert(x_buf.size() == static_cast<size_t>(order));
        assert(y_buf.size() == static_cast<size_t>(order));

        return true;
    };

    void verify_sample_rate(double dt) {
    	double dt_ = 1.0 / sample_rate;
        sample_rate_verified = std::fabs(dt - dt_) < 0.49 * dt_;
        if (!sample_rate_verified) {
        	char s[1024] = {0};
        	snprintf(s, 1024, "Sample rate does not match! desired interval[%f] input interval[%f]", 1.0 / sample_rate, dt);
            throw std::runtime_error(s);
        }
    };

  private:
    std::deque<double> x_buf;
    std::deque<double> y_buf;
    bool sample_rate_verified;
};

/* ==== MATLAB Script for filter design ====

fc = 5;
fs = 100;
order = 2;
[b,a] = butter(order,fc/(fs/2));

figure(1);
freqs(b,a);
figure(2);
grpdelay(b,a,100,fs);

*/

#include <iostream>
#include <random>
#include "matplot/matplot.h"
#include "nstd.h" 

#define EXERCISE 4
#define nstd_print(var) nstd::print(var, #var)

namespace ct
{
	constexpr double epsilon0 = 8.854187817e-12;
	constexpr double PI = 3.14159265359;
	constexpr double h = 6.62607015e-34;
	constexpr double hbar = ct::h / (2 * ct::PI);
}

struct charge
{
	// properties
private:
	double m_value = 1;
	double m_R = 1;
	// angles are in rad
	double m_phi;
	double m_theta;
	std::vector<double> m_position;

	// methods
public:
	charge(double phi, double theta, double R = 1, double value = 1) 
	{
		m_position.resize(3);
		updatePosition(phi, theta, R, value);
	}

	void updatePosition(double newPhi, double newTheta, double newR = 1, double newValue = 1)
	{
		m_phi = newPhi;
		m_theta = newTheta;
		m_R = newR;
		m_value = newValue;
		updateXYZPos();
	}

	double theta() const { return m_theta; }
	double phi() const { return m_phi; }
	double R() const { return m_R; }
	double value() const { return m_value; }
	const std::vector<double>& position() const { return m_position; }

private:
	void updateXYZPos()
	{
		m_position[0] = m_R * std::cos(m_phi) * std::sin(m_theta);
		m_position[1] = m_R * std::sin(m_phi) * std::sin(m_theta);
		m_position[2] = m_R * std::cos(m_theta);
	}
};

struct chargeGradient
{
	double m_R_gradient;
	double m_phi_gradient;
	double m_theta_gradient;
};

struct chargeDistribution
{
	// properties
public:
	double m_gradientStepSize = 1e-4;
	double m_gradientDescentIterations = 100000;

private:
	double m_energy;
	size_t m_numCharges;
	std::vector<charge> m_charges;
	std::vector<chargeGradient> m_gradient;

	// methods
public:
	chargeDistribution(std::vector<charge>& charges) 
	{
		m_numCharges = charges.size();
		m_gradient.resize(m_numCharges);
		updateDistribution(charges);
		calculatePotentialEnergy();
		findGradient();
	}

	void updateDistribution(std::vector<charge>& charges)
	{
		m_charges = charges;
		calculatePotentialEnergy();
	}

	void calculatePotentialEnergy()
	{
		m_energy = 0;
		double correctionFactor = 1;
		size_t numCharges = m_charges.size();
		std::vector<double> r_jk(numCharges);
		double r = 0;

		for (size_t i = 0; i < numCharges; ++i)
		{
			for (size_t j = (i + 1); j < numCharges; ++j)
			{
				charge& q1 = m_charges[i];
				charge& q2 = m_charges[j];

				double x1 = std::cos(q1.phi()) * std::sin(q1.theta());
				double y1 = std::sin(q1.phi()) * std::sin(q1.theta());
				double z1 = std::cos(q1.theta());

				double x2 = std::cos(q2.phi()) * std::sin(q2.theta());
				double y2 = std::sin(q2.phi()) * std::sin(q2.theta());
				double z2 = std::cos(q2.theta());

				//r = std::acos(std::sin(m_charges[i].theta()) * std::sin(m_charges[j].theta()) * std::cos(m_charges[i].phi() -
					//m_charges[j].phi()) + std::cos(m_charges[i].theta()) * std::cos(m_charges[j].theta()));
				r = std::sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
				m_energy += 1 / r * q1.value() * q2.value(); // * 1 / (4 * ct::PI * ct::epsilon0 * r_jk[i]);
			}
		}
	}

	void minimizeEnergy()
	{
		for (size_t i = 0; i < m_gradientDescentIterations; ++i)
		{
			findGradient();
			for (size_t j = 0; j < m_numCharges; ++j)
			{
				charge& q = m_charges[j];
				double theta = q.theta();
				double phi = q.phi();
				double R = q.R();
				theta += (-m_gradientStepSize * m_gradient[j].m_theta_gradient);
				phi += (-m_gradientStepSize * m_gradient[j].m_phi_gradient);

				q.updatePosition(phi, theta, R);
			}
		}
	}

	double energy() const { return m_energy; }
	std::vector<charge> charges() const { return m_charges; }
	std::vector<chargeGradient> gradient() const { return m_gradient; }

private:
	void findGradient()
	{
		double theta;
		double initialTheta;
		double phi;
		double initialPhi;
		double R;
		double initialR;
		double initialEnergy;

		for (size_t i = 0; i < m_charges.size(); ++i)
		{
			charge& q = m_charges[i];
			initialEnergy = m_energy;
			theta = q.theta();
			phi = q.phi();
			R = q.R();

			initialTheta = theta;
			initialPhi = phi;
			initialR = R;

			// theta
			theta += m_gradientStepSize;
			q.updatePosition(phi, theta, R);
			updateDistribution(m_charges);
			m_gradient[i].m_theta_gradient = (m_energy - initialEnergy) / m_gradientStepSize;
			theta = initialTheta;

			// phi
			phi += m_gradientStepSize;
			q.updatePosition(phi, theta, R);
			updateDistribution(m_charges);
			m_gradient[i].m_phi_gradient = (m_energy - initialEnergy) / m_gradientStepSize;
			phi = initialPhi;

			// R
			R += m_gradientStepSize;
			q.updatePosition(phi, theta, R);
			updateDistribution(m_charges);
			m_gradient[i].m_R_gradient = (m_energy - initialEnergy) / m_gradientStepSize;
			R = initialR;
			q.updatePosition(phi, theta, R);

		}

			updateDistribution(m_charges);
	}

};

static double potential(double x, double m, double omega)
{
	return 0.5 * m * m * omega * omega * x * x;
}

static double solveODE(double y_0, double dy_0, double x0, double L, double E, double dx, double m, double omega);

static std::vector<std::vector<double>> bracketing(double x0, double dx, double E, double L, double m, double omega)
{
	std::vector<std::vector<double>> brackets;
	double y0 = 0;
	double y1 = solveODE(y0, 1, x0, L, E, dx, m, omega);
	while (y1 * y0 > 0)
	{
		x0 += dx;
		y0 = y1;
		y1 = solveODE(y0, 1, x0, L, E, dx, m, omega);
	}
	brackets.push_back({ x0 - dx, x0 });
	brackets.push_back({ x0, x0 + dx });
	return brackets;
}

static std::vector<double> calcDerivatives(const std::vector<double>& y, double x, double E, double m, double omega)
{
	std::vector<double> derivatives(2);
	derivatives[0] = y[1];
	derivatives[1] = (2 * m / (ct::hbar * ct::hbar)) * (potential(x, m, omega) - E) * y[0];

	return derivatives;
}

static double solveODE(double y_0, double dy_0, double x0, double L, double E, double dx, double m, double omega)
{
	std::vector<double> y = { y_0, dy_0 };
	double x = x0;

	while (x <= L)
	{
		std::vector<double> dPsi = calcDerivatives(y, x, E, m, omega);

		y[0] += dx * dPsi[0];
		y[1] += dx * dPsi[1];

		x += dx;
	}

	return y[0];
}

static std::vector<double> findEnergies(double E, double L, double m, double omega, double tolerance, bool even = true)
{
	double x0 = 0;
	double y0 = 1;
	double dy0 = 0;
	if (!even) { y0 = 0; dy0 = 1; }
    std::vector<std::vector<double>> brackets = bracketing(x0, tolerance, E, L, m, omega);
    std::vector<double> energies;

    for (const auto& bracket : brackets)
    {
        y0 = solveODE(y0, dy0, bracket[0], L, E, tolerance, m, omega);
        double y1 = solveODE(y0, dy0, bracket[1], L, E, tolerance, m, omega);
        double x0 = bracket[0];
        double x1 = bracket[1];

        while (std::abs(x1 - x0) > tolerance)
        {
            double x2 = (x0 + x1) / 2;
            double y2 = solveODE(0, 1, x2, L, E, tolerance, m, omega);
            if (y0 * y2 < 0)
            {
                x1 = x2;
                y1 = y2;
            }
            else
            {
                x0 = x2;
                y0 = y2;
            }
        }
        energies.push_back((x0 + x1) / 2);
    }

    return energies;
}

int main()
{
#if EXERCISE == 1
	{
		const int numCharges = 4;
		std::vector<charge> charges;
		std::random_device rd;
		std::uniform_real_distribution<double> phiDist(0.0, 2 * ct::PI);
		std::uniform_real_distribution<double> thetaDist(0.0, ct::PI);

		for (int i = 0; i < numCharges; ++i)
		{
			charge newRandomCharge(phiDist(rd), thetaDist(rd));
			charges.push_back(newRandomCharge);
		}

		chargeDistribution dist(charges);
		dist.minimizeEnergy();

		std::vector<double> x1;
		std::vector<double> y1;
		std::vector<double> z1;

		for (size_t i = 0; i < dist.charges().size(); ++i)
		{
			x1.push_back(dist.charges()[i].position()[0]);
			y1.push_back(dist.charges()[i].position()[1]);
			z1.push_back(dist.charges()[i].position()[2]);
		}




#ifdef NDEBUG
		{
			using namespace matplot;
			const double R = 1.0;
			const int n = 50;

			std::vector<double> theta, phi;

			for (int i = 0; i < n; ++i) 
			{
				theta.push_back(ct::PI * i / (n - 1));
				phi.push_back(2 * ct::PI * i / (n - 1));
			}

			std::vector<std::vector<double>> X, Y, Z;

			for (double t : theta) 
			{
				std::vector<double> x_row, y_row, z_row;
				for (double p : phi) 
				{
					double x = R * sin(t) * cos(p);
					double y = R * sin(t) * sin(p);
					double z = R * cos(t);

					x_row.push_back(x);
					y_row.push_back(y);
					z_row.push_back(z);
				}
				X.push_back(x_row);
				Y.push_back(y_row);
				Z.push_back(z_row);
    }
		
			auto S = mesh(X, Y, Z);
			S->edge_color("k");
			S->hidden_3d(false);
			hold(on);
			plot3(x1, y1, z1, "r.");
			show();
		}
#endif
	}
#endif

#ifdef EXERCISE == 4
	{
		double m = 1.0, omega = 1.0;

		std::vector<double> k = findEnergies(0, 10, m, omega, 1e-6);
		nstd_print(k);

#ifdef NDEBUG
		{

		}
#endif
	}
#endif
}

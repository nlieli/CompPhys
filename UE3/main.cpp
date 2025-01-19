#include "matplot/matplot.h"
#include "nstd.h"
#include <iostream>
#include <random>
#include "SF.h"

#define EXERCISE 1
#define nstd_print(var) nstd::print(var, #var)

namespace ct
{
	constexpr double epsilon0 = 8.854187817e-12;
	constexpr double PI = 3.14159265359;
	constexpr double h = 6.62607015e-34;
	//constexpr double hbar = ct::h / (2 * ct::PI);
	constexpr double hbar = 1;
}

struct charge
{
	// properties
private:
	double m_value = 1;
	double m_R;
	// angles are in rad
	double m_phi;
	double m_theta;
	std::vector<double> m_position;

	// methods
public:
	charge(double phi, double theta, double R = 1, double value = 1) 
	{
		m_position.resize(3);
		updatePosition(phi, theta, R);
	}

	void updatePosition(double newPhi, double newTheta, double newR = 1)
	{
		m_phi = newPhi;
		m_theta = newTheta;
		m_R = newR;
		//m_value = newValue;
		updateXYZPos();
	}

	void changeChargeValue(double newValue)
	{
		m_value = newValue;
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
				m_energy += 1 / r * q1.value() * q2.value(); // *1 / (4 * ct::PI * ct::epsilon0);
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

struct mass
{
	double m_mass = 1;
	std::vector<double> position;

	mass() : position(3, 0.0) {}
};

struct t_mass : mass
{
	std::vector<double> velocity;
	std::vector<double> acceleration;

	t_mass() : mass(), velocity(3, 0.0), acceleration(3, 0.0) {}
};

struct C3BP
{
	double dt = 1e-3;

	mass m_M1;
	mass m_M2;
	t_mass m_m;
	
	std::vector<double> omega = { 0, 0, 1 };
	double mu;

	std::vector<double> r1;
	std::vector<double> r2;

	std::vector<double> F_g;
	std::vector<double> F_i;
	std::vector<double> F_t;

	std::vector<std::vector<double>> trajectory;

	C3BP(const mass& M1, const mass& M2, const t_mass& m) : m_M1(M1), m_M2(M2), m_m(m)
	{
		mu = m_M2.m_mass / (m_M1.m_mass + m_M2.m_mass);
		m_M1.position[0] = -mu;
		m_M2.position[0] = 1 - mu;
		trajectory.resize(3);

		//trajectory[0].push_back(m_m.position[0]);
		//trajectory[1].push_back(m_m.position[1]);
		//trajectory[2].push_back(m_m.position[2]);

		//calculateForces();
	}

	void calculateForces()
	{
		r1 = m_m.position - m_M1.position;
		r2 = m_m.position - m_M2.position;

		F_g = (mu - 1) * r1 / (std::pow(nstd::norm1dimArray(r1), 3)) - mu * r2 / (std::pow(nstd::norm1dimArray(r2), 3));
		std::vector<double> temp_vec = nstd::crossProduct(omega, m_m.position);
		F_i = -2 * nstd::crossProduct(omega, m_m.velocity) - nstd::crossProduct(omega, temp_vec);
		F_t = F_g + F_i;
	};

	inline void updatePosition()
	{
		m_m.acceleration = F_t;
		m_m.velocity = m_m.velocity + m_m.acceleration * dt;
		m_m.position = m_m.position + m_m.velocity * dt;
		trajectory[0].push_back(m_m.position[0]);
		trajectory[1].push_back(m_m.position[1]);
		trajectory[2].push_back(m_m.position[2]);
		calculateForces();
	}

	inline void updatePositionRK4()
	{
		auto position = m_m.position;
		auto velocity = m_m.velocity;

		calculateForces();
		auto k1_v = F_t;                       
		auto k1_p = velocity;                  

		auto temp_position = position + k1_p * (dt / 2.0);
		auto temp_velocity = velocity + k1_v * (dt / 2.0);
		m_m.position = temp_position;
		m_m.velocity = temp_velocity;
		calculateForces();
		auto k2_v = F_t;
		auto k2_p = temp_velocity;

		temp_position = position + k2_p * (dt / 2.0);
		temp_velocity = velocity + k2_v * (dt / 2.0);
		m_m.position = temp_position;
		m_m.velocity = temp_velocity;
		calculateForces();
		auto k3_v = F_t;
		auto k3_p = temp_velocity;

		temp_position = position + k3_p * dt;
		temp_velocity = velocity + k3_v * dt;
		m_m.position = temp_position;
		m_m.velocity = temp_velocity;
		calculateForces();
		auto k4_v = F_t;
		auto k4_p = temp_velocity;

		m_m.position = position + (k1_p + k2_p * 2.0 + k3_p * 2.0 + k4_p) * (dt / 6.0);
		m_m.velocity = velocity + (k1_v + k2_v * 2.0 + k3_v * 2.0 + k4_v) * (dt / 6.0);

		trajectory[0].push_back(m_m.position[0]);
		trajectory[1].push_back(m_m.position[1]);
		trajectory[2].push_back(m_m.position[2]);
	}

	void createTrajectory(int numberOfTimeSteps)
	{
		for (int i = 0; i < numberOfTimeSteps; ++i)
			updatePositionRK4();
	}

	std::vector<std::vector<double>> lagrangePoints() const
	{
		std::vector<std::vector<double>> LP(3, std::vector<double>(5, 0.0));
		std::vector<std::vector<double>> quintics = {
			{ -mu, 2 * mu, -mu, 3 - 2 * mu, mu - 3, 1 },
			{ -mu, -2 * mu, -mu, 3 - 2 * mu, 3 - mu, 1},
			{ -mu, 12 + 14 * mu, -24 - 13 * mu, 6 * mu + 19, -7 - mu, 1} };

		for (size_t i = 0; i < quintics.size(); ++i)
		{
			std::vector<double> intervals = polyBracketing(quintics[i], -5, 5);
			LP[0][i] = polyNewtonRaphson(quintics[i], 10, intervals[0]);
		}

		LP[0][0] = (1 - mu) - LP[0][0];
		LP[0][1] += 1 - mu;
		LP[0][2] = -1 - LP[0][2];

		LP[0][3] = std::cos(ct::PI / 3) - mu;
		LP[1][3] = std::sin(ct::PI / 3);

		LP[0][4] = std::cos(-ct::PI / 3) - mu;
		LP[1][4] = std::sin(-ct::PI / 3);

		return LP;
	}

	void setTestinLagrangePoint_x(int x)
	{
		if (x < 1 || x > 5) { std::cerr << "\033[1;31m[ERROR]: Lagrange Point must be int from 1 to 5.\033[0m"; throw std::runtime_error("Not a lagrange point"); }
		--x;
		std::vector<std::vector<double>> L = lagrangePoints();
		m_m.position = { L[0][x], L[1][x], L[2][x] };
		nstd_print(m_m.position);
	}
};

static double potential(double x, double m, double omega, double lambda)
{
	return 0.5 * m * m * omega * omega * x * x + lambda / 24 * std::pow(x, 4);
}

static double solveODE(double y_0, double dy_0, double x0, double L, double E, double m, double omega, double lambda);

static std::vector<std::vector<double>> bracketing(double y_0, double dy_0, double x0, double L, double m, double omega, double lambda, int N)
{
	std::vector<std::vector<double>> brackets;
	double E0 = 0;
	double E1;
	double dE = 1e-1;

	double psi0 = solveODE(y_0, dy_0, x0, L, E0, m, omega, lambda);
	double psi1;

	int i = 0;

	while (i < N)
	{
		E1 = E0 + dE;
		psi1 = solveODE(y_0, dy_0, x0, L, E1, m, omega, lambda);

		if (psi1 * psi0 < 0)
		{
			brackets.push_back({ E0, E1 });
			++i;
		}
		psi0 = psi1;
		E0 = E1;
	}

	return brackets;
}

static std::vector<double> calcDerivatives(std::vector<double> y, double x, double E, double m, double omega, double lambda)
{
	std::vector<double> derivatives(2);
	derivatives[0] = y[1];
	derivatives[1] = 2 * m * (potential(x, m, omega, lambda) - E) * y[0] / (ct::hbar * ct::hbar);

	return derivatives;
}

static double solveODE(double y_0, double dy_0, double x0, double L, double E, double m, double omega, double lambda)
{
	std::vector<double> y = { y_0, dy_0 };
	double x = x0;
	double dx = 1e-3;
	std::vector<double> dPsi;

	while (x <= L)
	{
		// using RK4
        std::vector<double> k1 = calcDerivatives(y, x, E, m, omega, lambda);

        std::vector<double> y_temp1 = { y[0] + 0.5 * dx * k1[0], y[1] + 0.5 * dx * k1[1] };
        std::vector<double> k2 = calcDerivatives(y_temp1, x + 0.5 * dx, E, m, omega, lambda);

        std::vector<double> y_temp2 = { y[0] + 0.5 * dx * k2[0], y[1] + 0.5 * dx * k2[1] };
        std::vector<double> k3 = calcDerivatives(y_temp2, x + 0.5 * dx, E, m, omega, lambda);

        std::vector<double> y_temp3 = { y[0] + dx * k3[0], y[1] + dx * k3[1] };
        std::vector<double> k4 = calcDerivatives(y_temp3, x + dx, E, m, omega, lambda);

        y[0] += (dx / 6.0) * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]);
        y[1] += (dx / 6.0) * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]);

        x += dx;
	}

	return y[0];
}

static std::vector<double> findEnergies(double y_0, double dy_0, double x_0, double L, double m, double omega, double lambda, double tolerance, double NumberOfValues)
{
	std::vector<std::vector<double>> brackets = bracketing(y_0, dy_0, x_0, L, m, omega, lambda, NumberOfValues);

	std::vector<double> energies;

	double El, Eu;
	double psiL, psiU;

	for (const std::vector<double>& bracket : brackets)
	{
		El = bracket[0];
		Eu = bracket[1];

		psiL = solveODE(y_0, dy_0, x_0, L, El, m, omega, lambda);
		psiU = solveODE(y_0, dy_0, x_0, L, Eu, m, omega, lambda);

		while (std::abs(Eu - El) > tolerance)
		{
			double Em = 0.5 * (El + Eu);
			double psiM = solveODE(y_0, dy_0, x_0, L, Em, m, omega, lambda);
			if (psiM * psiL < 0)
			{
				Eu = Em;
				psiU = psiM;
			}
			else
			{
				El = Em;
				psiL = psiM;
			}
		}

		energies.push_back(0.5 * (El + Eu));
	}

	return energies;
}

int main()
{
#if EXERCISE == 1
	{
		const int numCharges = 6;
		std::vector<charge> charges;
		std::random_device rd;
		std::uniform_real_distribution<double> phiDist(0.0, 2 * ct::PI);
		std::uniform_real_distribution<double> thetaDist(0.0, ct::PI);

		for (int i = 0; i < numCharges; ++i)
		{
			charge newRandomCharge(phiDist(rd), thetaDist(rd));
			charges.push_back(newRandomCharge);
		}

		//charges[0].changeChargeValue(100);
		
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

#if EXERCISE == 3
	{
		mass M1;
		// M1.position = { 0.5, 0, 0 }; // irrelevant
		M1.m_mass = 1e7;
		mass M2;
		// M2.position = { -0.5, 0, 0 }; // irrelevant
		M2.m_mass = 1e5;
		t_mass m;
		m.position = { 0.4, 0.8, 0 };
		m.m_mass = 1;
		m.velocity = { 0, 0.02, 0 };

		C3BP system(M1, M2, m);
		system.omega[2] = 1;
		system.dt = 1e-3;
		std::vector<std::vector<double>> L = system.lagrangePoints();
		system.setTestinLagrangePoint_x(4);

		system.createTrajectory(50000);

#ifdef NDEBUG
		{
			using namespace matplot;
			std::vector<double>& M1_pos = system.m_M1.position;
			std::vector<double>& M2_pos = system.m_M2.position;

			std::vector<double> M1x = { M1_pos[0] };
			std::vector<double> M1y = { M1_pos[1] };
			std::vector<double> M1z = { M1_pos[2] };

			std::vector<double> M2x = { M2_pos[0] };
			std::vector<double> M2y = { M2_pos[1] };
			std::vector<double> M2z = { M2_pos[2] };

			std::vector<double>& Tx = system.trajectory[0];
			std::vector<double>& Ty = system.trajectory[1];
			std::vector<double>& Tz = system.trajectory[2];

			std::vector<double> Ix = { system.trajectory[0][0] };
			std::vector<double> Iy = { system.trajectory[1][0] };
			std::vector<double> Iz = { system.trajectory[2][0] };

			figure();
			auto K = plot3(M1x, M1y, M1z, ".r");
			K->line_width(2.0);
			hold(on);
			plot3(M2x, M2y, M2z, ".g");
			plot3(L[0], L[1], L[2], ".k");
			plot3(Ix, Iy, Iz, "ob");
			plot3(Tx, Ty, Tz, "b");

			double lim = 3;
			xlim({ -lim, lim });
			ylim({ -lim, lim });
			zlim({ -1, 1 });
			auto S = legend("M1", "M2", "Lagrange Points");
			//S->position({ -5.0, -5.0 });
			show();
		}
#endif
	}
#endif

#if EXERCISE == 4
	{
		double m = 1.0, omega = 1.0;
		double L = 10, tolerance = 1e-10;
		double lambda = 24;
		auto E = findEnergies(0, 1, 0, L, m, omega, lambda, tolerance, 5);
		auto E1 = findEnergies(1, 0, 0, L, m, omega, lambda, tolerance, 5);
		nstd_print(E);
		nstd_print(E1);


#ifdef NDEBUG
		{

		}
#endif
	}
#endif
}

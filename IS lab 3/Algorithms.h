#pragma once
#include <random>
#include <iostream>
#include <unordered_set>
#include <set>
#include <bitset>
#include <sstream>
const double PI = 3.141592653589793238463;

template <typename T>
struct Random {
	/// <summary>
	/// Возвращает случайное число от 0.0 до 1.0 включительно
	/// </summary>
	static T random() {
		return static_cast<T>(rand()) / static_cast<T>(RAND_MAX);
	}
	/// <summary>
	/// Возвращает случайное число от 'from' до 'to' включительно
	/// </summary>
	static T random(T from, T to) {
		return from + static_cast<T>(rand()) / (static_cast<T>(RAND_MAX) / static_cast<T>(to - from));
	}
};

template <typename T>
struct State {
	T first;
	T second;
	State(T first, T second) {
		this->first = first;
		this->second = second;
	}
	T& operator[](size_t ind) {
		if (ind == 0) {
			return first;
		}
		return second;
	}
	const T& operator[](size_t ind) const {
		if (ind == 0) {
			return first;
		}
		return second;
	}
	/// <summary>
	/// Возвращает случайное состояние в заданных пределах
	/// </summary>
	static State<T> gen_random_state(T from, T to) {
		T first = Random<int>::random(from, to);
		T second = Random<int>::random(from, to);
		return State<T>(first, second);
	}
};

enum ExtremumType{min, max};

template <typename T>
class Annealing
{
	T(*func)(T, T);
	T temp_coef = 0.97;
	T start_temp;
	T end_temp;
	T step_size;
	
	/// <summary>
	/// Возвращает энергию системы в состоянии
	/// </summary>
	T E(State<T> state) {
		return (*func)(state.first, state.second);
	}
	/// <summary>
	/// Возвращает вероятность принятия нового состояния
	/// </summary>
	T p(T delta_e, T temperature) {
		return std::exp(-delta_e / temperature);
	}
	/// <summary>
	/// Возвращает надо ли принимать новое состояние
	/// </summary>
	bool accept_state(T delta_e, T temperature) {
		T border = p(delta_e, temperature);
		T alpha = Random<T>::random();
		return alpha <= border;
	}
	/// <summary>
	/// Возвращает температуру на k-ом шаге по схеме больцмановского отжига
	/// </summary>
	T boltzmann_temperature(unsigned int k){
		return start_temp / std::log(1 + k);
	}
	/// <summary>
	/// Возвращает температуру на k-ом шаге по схеме тушения
	/// </summary>
	T extinguishing_temperature(T last_temp) {
		return last_temp * temp_coef;
	}
	/// <summary>
	/// Возвращает Евклидово расстояние между двумя состояниями
	/// </summary>
	T euclidean_dist(State<T> first_state, State<T> second_state) {
		return std::sqrt(
			std::pow(first_state.first - second_state.first, 2) +
			std::pow(first_state.second - second_state.second, 2)
			);
	}
	/// <summary>
	/// Порождающее семейство вероятностных распределений ζ(x,T)
	///  выбирается как семейство нормальных распределений 
	/// с математическим ожиданием x (old_state) и дисперсией T (temperature)
	/// </summary>
	T g(State<T> old_state, T temperature) {
		// изначально первой степенью была -n/2, где n - размерность функции, так что пусть остается
		T part = 1.0 / 2 * PI * temperature;
		//T elem = part *
		//	std::exp(-std::pow(-std::abs(old_state.second - old_state.first), 2) / (2 * temperature));
		//elem = part;
		return part;
	}
	/// <summary>
	/// Порождающее семейство вероятностных распределений ζ(x,T)
	///  выбирается как семейство нормальных распределений 
	/// с математическим ожиданием x (old_state) и дисперсией T (temperature)
	/// </summary>
	State<T> gen_new_state(State<T> old_state, T temperature) {
		auto step = g(old_state, temperature);
		int coef = 1;
		if (Random<T>::random() < 0.5) {
			coef = -1;
		}
		T first = old_state.first + step * coef;
		coef = 1;
		if (Random<T>::random() < 0.5) {
			coef = -1;
		}
		T second = old_state.second + step * coef;
		return State<T>(first, second);
	}
public:
	Annealing(T start_temp, T end_temp, T step_size, T(*func)(T, T)) {
		this->start_temp = start_temp;
		this->end_temp = end_temp;
		this->step_size = step_size;
		this->func = func;
	}
	State<T> find(T from, T to, ExtremumType ex_type) {
		auto state = State<T>::gen_random_state(from, to);	// случайное начальное состояние
		T e = E(state);							// начальная энергия
		unsigned int k = 1;					// текущий шаг
		T tk = start_temp;					// максимальная начльная температура
		int ex_multiplier;
		if (ex_type == min) {
			ex_multiplier = 1;
		}
		else {
			ex_multiplier = -1;
		}
		while (tk >= end_temp) {
			while (true) {
				auto new_state = gen_new_state(state, tk);
				T new_e = E(new_state);
				T delta_e = (new_e - e) * ex_multiplier;
				if (delta_e >= 0 && !accept_state(delta_e, tk)) {
					continue;
				}
				state = new_state;
				e = new_e;
				k += 1;
				tk = extinguishing_temperature(tk);
				break;
			}
		}

		return state;
	}
};

template <typename T>
class Genetic
{
	T(*func)(T, T);
	size_t population_size;
	size_t max_epoch;

	/// <summary>
	/// Возвращает набор начальных состояний
	/// </summary>
	std::vector<State<T>> gen_start_population(T from, T to) {
		std::vector<State<T>> pops;
		for (size_t i = 0; i < population_size; i++) {
			pops.push_back(State<T>::gen_random_state(from, to));
		}
		return pops;
	}
	/// <summary>
	/// Возвращает набор значений функции для переданных состояний
	/// </summary>
	std::vector<T> estimate(const std::vector<State<T>>& pops) {
		std::vector<T> vals;
		for (const State<T>& pop : pops) {
			vals.push_back((*func)(pop.first, pop.second));
		}
		return vals;
	}
	/// <summary>
	/// Возвращает индекс минимального элемента
	/// </summary>
	size_t get_ind_with_ex_val(const std::vector<T>& vals, int ex_multiplier) {
		size_t ex_ind = 0;
		T ex = vals[0];
		for (size_t i = 1; i < vals.size(); i++) {
			if (vals[i] * ex_multiplier < ex * ex_multiplier) {
				ex = vals[i];
				ex_ind = i;
			}
		}
		return ex_ind;
	}
	/// <summary>
	/// Возвращает знаки элментов состояний
	/// </summary>
	std::vector<int> get_signs(const std::vector<State<T>>& states) {
		std::vector<int> signs;
		for (size_t i = 0; i < 2; i++) {
			int sign;
			if (states[0][i] < 0 && states[1][i] < 0) {
				sign = -1;
			}
			else if (states[0][i] >= 0 && states[1][i] >= 0) {
				sign = 1;
			}
			else {
				sign = Random<int>::random(-1, 1);
			}
			signs.push_back(sign);
		}

		return signs;
	}
	/// <summary>
	/// Возвращает состояния, закодированные в битовое представление 
	/// </summary>
	std::vector<State<std::bitset<64>>> encode_states(const std::vector<State<T>>& states) {
		std::vector<State<std::bitset<64>>> bits_states;
		for (size_t i = 0; i < 2; i++) {
			auto first = std::bitset<64>(*reinterpret_cast<const uint64_t*>(&states[i].first));
			auto second = std::bitset<64>(*reinterpret_cast<const uint64_t*>(&states[i].second));
			bits_states.push_back(State<std::bitset<64>>(first, second));
		}
		return bits_states;
	}

	union Converter { uint64_t i; double d; };
	double bits_to_double(std::bitset<64> const& bs) {
		Converter c;
		c.i = bs.to_ullong();
		return c.d;
	}
	double bitstring_to_double(const std::string& s)
	{
		unsigned long long x = 0;
		for (std::string::const_iterator it = s.begin(); it != s.end(); ++it)
		{
			x = (x << 1) + (*it - '0');
		}
		double d;
		memcpy(&d, &x, 8);
		return d;
	}
	/// <summary>
	/// Возвращает состояние из закодированного в битовое представление состояния
	/// </summary>
	State<T> decode_state(const State<std::bitset<64>>& encoded_states,
		const std::vector<int>& signs) {
		auto decoded_state = State<T>(bits_to_double(encoded_states[0]),
			bits_to_double(encoded_states[1]));
		return decoded_state;
	}
	/// <summary>
	/// Возвращает два новых закодированных состояния на основе переданных
	/// </summary>
	State<std::bitset<64>> cross_states(const std::vector<State<std::bitset<64>>>& encoded_states) {
		auto elem = std::bitset<64>(0);
		auto crossed_state = State<std::bitset<64>>(elem, elem);

		auto state_size = encoded_states[0][0].size();
		auto mutation_ver = 1.0 / state_size;
		for (size_t i = 0; i < 2; i++) {
			auto fes = encoded_states[0][i];
			auto ses = encoded_states[1][i];
			for (size_t j = 0; j < state_size; j++) {
				if (fes[j] == ses[j]) {
					if (Random<double>::random() <= mutation_ver) {
						crossed_state[i].set(j, ~fes[j]);
						continue;
					}
					crossed_state[i].set(j, fes[j]);
				}
				else {
					auto rnd = Random<int>::random(0, 1);
					crossed_state[i].set(j, rnd);
				}
			}
		}

		return crossed_state;
	}
	void next_generation(std::vector<State<T>>& pops, int ex_multiplier) {
		auto vals = estimate(pops);
		auto ind = get_ind_with_ex_val(vals, ex_multiplier);
		std::vector<State<T>> children;
		children.push_back(pops[ind]);

		while (children.size() < population_size) {
			std::vector<State<T>> states;
			for (size_t i = 0; i < 2; i++) {
				auto f_ind = Random<size_t>::random(0, population_size - 1);
				auto s_ind = Random<size_t>::random(0, population_size - 1);
				if (vals[f_ind] * ex_multiplier < vals[s_ind] * ex_multiplier) {
					states.push_back(pops[f_ind]);
				}
				else {
					states.push_back(pops[s_ind]);
				}
			}

			auto signs = get_signs(states);
			auto bits_states = encode_states(states);
			auto crossed_state = cross_states(bits_states);
			auto decoded_state = decode_state(crossed_state, signs);

			children.push_back(decoded_state);
		}

		pops = children;
	}

public:
	Genetic(size_t population_size, size_t max_epoch, T(*func)(T, T)) {
		this->population_size = population_size;
		this->max_epoch = max_epoch;
		this->func = func;
	}
	State<T> find(T from, T to, ExtremumType ex_type) {
		auto pops = gen_start_population(from, to);
		int ex_multiplier;
		if (ex_type == min) {
			ex_multiplier = 1;
		}
		else {
			ex_multiplier = -1;
		}
		for (size_t i = 0; i < max_epoch; i++) {
			next_generation(pops, ex_multiplier);
		}
		return pops[0];
	}
};

template <typename T>
class HillClimbing {
	T(*func)(T, T);
	T step_size;

	/// <summary>
	/// Возвращает возможные новые состояния
	/// </summary>
	std::vector<State<T>> gen_candidates(const State<T>& state) {
		std::vector<State<T>> candidates;
		candidates.push_back(State<T>(state.first - step_size, state.second));
		candidates.push_back(State<T>(state.first - step_size, state.second - step_size));
		candidates.push_back(State<T>(state.first - step_size, state.second + step_size));
		candidates.push_back(State<T>(state.first + step_size, state.second));
		candidates.push_back(State<T>(state.first + step_size, state.second - step_size));
		candidates.push_back(State<T>(state.first + step_size, state.second + step_size));
		candidates.push_back(State<T>(state.first, state.second - step_size));
		candidates.push_back(State<T>(state.first, state.second + step_size));
		return candidates;
	}
	/// <summary>
	/// Возвращает новое состояние в случайном направлении
	/// </summary>
	State<T> gen_next_state(const State<T>& state) {
		int coef = 1;
		if (Random<T>::random() < 0.5) {
			coef = -1;
		}
		T first = state.first + step_size * coef;
		coef = 1;
		if (Random<T>::random() < 0.5) {
			coef = -1;
		}
		T second = state.second + step_size * coef;
		return State<T>(first, second);
	}
	/// <summary>
	/// Возвращает новое состояние по жадному алгоритму
	/// </summary>
	bool greedy_set_state(State<T>& state, const std::vector<State<T>>& candidates, double ex_multiplier) {
		State<T> best_state = state;
		T val = (*func)(state.first, state.second);
		T best_val = val;
		for (size_t i = 0; i < candidates.size(); i++) {
			T cand_val = (*func)(candidates[i].first, candidates[i].second);
			if (cand_val * ex_multiplier < best_val * ex_multiplier) {
				best_val = cand_val;
				best_state = candidates[i];
			}
		}
		state = best_state;
		return (std::abs(val - best_val) < 1e-15);
	}
	/// <summary>
	/// Возвращает новое состояние по cтохастический алгоритму
	/// </summary>
	bool stochastic_set_state(State<T>& state, const std::vector<State<T>>& candidates, double ex_multiplier) {
		State<T> best_state = state;
		T val = (*func)(state.first, state.second);
		T best_val = val;

		std::unordered_set<size_t> set;
		while (set.size() != candidates.size()) {
			auto i = Random<size_t>::random(0, candidates.size() - 1);
			if (set.find(i) == set.end()) {
				T cand_val = (*func)(candidates[i].first, candidates[i].second);
				if (cand_val * ex_multiplier < best_val * ex_multiplier) {
					best_val = cand_val;
					best_state = candidates[i];
					break;
				}
			}
			set.insert(i);
		}
		state = best_state;
		return (std::abs(val - best_val) < 1e-15);
	}
	/// <summary>
	/// Возвращает новое состояние с выбором первого варианта
	/// </summary>
	bool first_set_state(State<T>& state, const std::vector<State<T>>& candidates, double ex_multiplier) {
		State<T> best_state = state;
		T val = (*func)(state.first, state.second);
		T best_val = val;
		for (size_t i = 0; i < candidates.size(); i++) {
			T cand_val = (*func)(candidates[i].first, candidates[i].second);
			if (cand_val * ex_multiplier < best_val * ex_multiplier) {
				best_val = cand_val;
				best_state = candidates[i];
				break;
			}
		}
		state = best_state;
		return (std::abs(val - best_val) < 1e-15);
	}
public:
	HillClimbing(T step_size, T(*func)(T, T)) {
		this->step_size = step_size;
		this->func = func;
	}
	State<T> find(T from, T to, ExtremumType ex_type) {
		auto state = State<T>::gen_random_state(from, to);

		double ex_multiplier;
		if (ex_type == min) {
			ex_multiplier = 1;
		}
		else {
			ex_multiplier = -1;
		}

		bool end = false;
		while (!end) {
			auto candidates = gen_candidates(state);
			end = greedy_set_state(state, candidates, ex_multiplier);
		}
		return state;
	}
	State<T> find_stochastic(T from, T to, ExtremumType ex_type) {
		auto state = State<T>::gen_random_state(from, to);

		double ex_multiplier;
		if (ex_type == min) {
			ex_multiplier = 1;
		}
		else {
			ex_multiplier = -1;
		}

		bool end = false;
		while (!end) {
			auto candidates = gen_candidates(state);
			end = stochastic_set_state(state, candidates, ex_multiplier);
		}
		return state;
	}
	State<T> find_first(T from, T to, ExtremumType ex_type) {
		auto state = State<T>::gen_random_state(from, to);

		double ex_multiplier;
		if (ex_type == min) {
			ex_multiplier = 1;
		}
		else {
			ex_multiplier = -1;
		}

		bool end = false;
		while (!end) {
			auto candidates = gen_candidates(state);
			end = first_set_state(state, candidates, ex_multiplier);
		}
		return state;
	}
	State<T> find_first_with_restart(T from, T to, ExtremumType ex_type) {
		double ex_multiplier;
		if (ex_type == min) {
			ex_multiplier = 1;
		}
		else {
			ex_multiplier = -1;
		}

		auto best_state = find_first(from, to, ex_type);
		auto best_val = (*func)(best_state.first, best_state.second);
		for (size_t i = 1; i < 500; i++) {
			auto state = find_first(from, to, ex_type);
			auto state_val = (*func)(state.first, state.second);
			T delta = state_val - best_val;
			if (delta * ex_multiplier < 0) {
				best_state = state;
				best_val = state_val;
			}
		}
		return best_state;
	}
	State<T> find_ray(T from, T to, ExtremumType ex_type) {
		auto state = State<T>::gen_random_state(from, to);

		double ex_multiplier;
		if (ex_type == min) {
			ex_multiplier = 1;
		}
		else {
			ex_multiplier = -1;
		}

		std::vector<State<T>> front = gen_candidates(state);
		State<T> best = front[0];
		bool end = false;
		while (!end) {
			std::vector<State<T>> front_candidates;
			for (size_t i = 0; i < front.size(); i++) {
				auto candidates = gen_candidates(front[i]);
				front_candidates.insert(front_candidates.end(), candidates.begin(), candidates.end());
			}
			auto func = this->func;
			std::sort(front_candidates.begin(), front_candidates.end(),
				[func, ex_multiplier](State<T> const& a, State<T> const& b) { 
					return (*func)(a[0], a[1]) * ex_multiplier < (*func)(b[0], b[1])* ex_multiplier; }
			);
			front = std::vector<State<T>>(front_candidates.begin(), front_candidates.begin() + front.size());
			if (std::abs(front[0].first - best.first) < step_size && 
				std::abs(front[0].second - best.second) < step_size) {
				break;
			}
			best = front[0];
		}
		return best;
	}
};
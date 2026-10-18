// Copyright @ Chun Shen 2018

#ifndef SRC_GLAUBER_H_
#define SRC_GLAUBER_H_

#include "data_structs.h"
#include "Parameters.h"
#include "RandomUlty.h"
#include "Nucleus.h"
#include "CollisionEvent.h"
#include "QCDString.h"
#include <memory>
#include <string>
#include <set>
#include <vector>

using std::shared_ptr;
class MCGlauberWrapper;
namespace MCGlb {

class Glauber {
 private:
    const Parameters &parameter_list;
    std::unique_ptr<Nucleus> projectile;
    std::unique_ptr<Nucleus> target;
    std::set<shared_ptr<CollisionEvent>, compare_collision_time> collision_schedule;

    std::vector<QCDString> QCD_string_list;
    std::vector<QCDString> remnant_string_list_;
    std::vector<CollisionEvent> collision_schedule_list_;
    std::vector<std::vector<real>> QCD_string_output_arr_;
    std::shared_ptr<RandomUtil::Random> ran_gen_ptr_;
    std::vector<double> Proj_nucleonz_;
    std::vector<double> Targ_nucleonz_;
    std::vector<double> HardPartonPosAndMomProj_;
    std::vector<double> HardPartonPosAndMomTarg_;
    std::vector<double> HardPartonPos_;
    std::vector<double> Mom_remnant_proj_;
    std::vector<double> Mom_remnant_targ_;
    bool sample_valence_quark;
    bool fluct_Nstrings_per_NN_collision_;
    real remnant_energy_loss_fraction_;

    std::vector<std::vector<real>> participantList_;

    real impact_b;
    real yloss_param_slope;
    real yloss_param_a;
    real yloss_param_b;

    real ybeam;

    int system_status_;

    real sigma_eff_;
    real nucleon_width_;

    SpatialVec Proj_hot_spot_x_;
    SpatialVec Targ_hot_spot_x_;
 public:
    //Glauber() = default;
    Glauber(const MCGlb::Parameters &param_in,
            shared_ptr<RandomUtil::Random> ran_gen);
    ~Glauber() {};

    std::vector<CollisionEvent> get_collision_information() {
        return (collision_schedule_list_);
    }

    void make_nuclei();
    real get_impact_parameter() const {return(impact_b);}

    int make_collision_schedule();
    int make_collision_schedule_second_stage();
    bool hit(real d2) const;

    int get_Npart() const;

    //! This function creates a new collision event between two nucleons
    void create_a_collision_event(shared_ptr<Nucleon> proj,
                                  shared_ptr<Nucleon> targ);
    bool get_collision_point(real t, real z1, real v1, real z2, real v2,
                             real &t_coll, real &z_coll) const;

    real compute_NN_inelastic_cross_section(real ecm) const;
    real GG_probality(real x, real Lambda, real Omega, real Sigma0) const;

    //! this function decides which of those binary collisions will produce
    //! QCD strings
    int decide_QCD_strings_production();
    int decide_QCD_strings_production_second_stage();
    //! this function determines whether a given binary collision event
    //! will produce a string
    int decide_produce_string_num(shared_ptr<CollisionEvent> event_ptr) const;

    void sample_valence_quarks();
    real get_tau_form(const int string_evolution_mode) const;
    void get_tau_form_and_moversigma(const int string_evolution_mode,
                                     const real y_in_lrf,
                                     real &tau_form, real &m_over_sigma,
                                     real &y_loss);

    real sample_rapidity_loss_shell(real y_init) const;
    real sample_rapidity_loss_from_the_LEXUS_model(const real y_init) const;
    real sample_rapidity_loss_from_parametrization(const real y_init) const;
    real sample_rapidity_loss_from_piecewise_parametrization(
                                                const real y_init) const;
    real sample_rapidity_loss_from_parametrization_with_fluct(
                                                const real y_init) const;

    real sample_junction_rapidity_right(real y_left, real y_right) const;
    real sample_junction_rapidity_left(real y_left, real y_right) const;
    real sample_junction_rapidity_uniformed(real y_left, real y_right) const;

    //! This function gets the target/projectile nucleon density at Lab frame
    //! at t, x, y, z The unit is 1/fm^3
    double get_nucleus_density(double t, double x, double y, double z,
                               int direction,
                               std::unique_ptr<Nucleus> &nucleus_ptr);

    double get_targ_nucleon_density(double t, double x, double y, double z);
    double get_proj_nucleon_density(double t, double x, double y, double z);

    //! This function gets the total nucleon density at Lab frame at t, x, y, z
    //! The unit is 1/fm^3
    double get_nucleon_density(double t, double x, double y, double z);

    std::vector<double> get_all_proj_nucleon_z();
    std::vector<double> get_all_targ_nucleon_z();

    std::vector<double>  OutputquarkPosProj();
    std::vector<double>  OutputquarkPosTarg();

    std::vector<double>  GetRemMomTarg() {
        return(Mom_remnant_targ_);
    }
    std::vector<double>  GetRemMomProj() {
        return(Mom_remnant_proj_);
    }

    void set_Proj_hot_spot_x(SpatialVec x) {Proj_hot_spot_x_ = x;}
    void set_Targ_hot_spot_x(SpatialVec x) {Targ_hot_spot_x_ = x;}

    SpatialVec get_Proj_hot_spot_x() const {return(Proj_hot_spot_x_);}
    SpatialVec get_Targ_hot_spot_x() const {return(Targ_hot_spot_x_);}

    //! This function performs string production between each nucleon pair
    int perform_string_production();
    void produce_remnant_strings();

    //! This function propagate individual nucleon inside the nucleus by dt
    void propagate_nuclei(real dt);
    void propagate_nucleon(shared_ptr<Nucleon> n_i, real dt);
    void update_momentum(shared_ptr<Nucleon> n_i, real y_shift);
    void update_momentum_quark(shared_ptr<Quark> q_i, real y_shift);
    //! This function updates the collision schedule
    void update_collision_schedule(shared_ptr<CollisionEvent> event_happened);

    void prepare_output_QCD_strings();
    //! the output rows of the strings and remnant strings, centred on (x_o, y_o)
    void fill_QCD_string_output(std::vector<QCDString> &strings,
                                std::vector<QCDString> &remnants,
                                const real x_o, const real y_o,
                                std::vector<std::vector<real>> &out) const;
    //! LO hard subtraction, mode 1: take E_take from the valence quarks of
    //! nucleon n (on the projectile side if proj_side) in strings, a copy of
    //! the event's; fills info's entries of side i_side (0 proj, 1 targ)
    void take_energy_from_quarks(std::vector<QCDString> &strings,
                                 shared_ptr<Nucleon> n, const bool proj_side,
                                 const real E_take, const int i_side,
                                 std::vector<real> &info) const;
    void computeCenterOfMass(real &x_o, real &y_o);
    void output_QCD_strings(std::string filename, const real Npart,
                            const real Ncoll, const real Nstrings,
                            const real b, const unsigned int seed);

    //void Pick_and_subtract_hard_parton_momentum_in_nucleon();
    void Pick_and_subtract_hard_parton_momentum();
    //! the target's hard quark: a position (relative to its nucleon at
    //! targ_n) from the valence-quark Gaussian (BG) within 0.5 fm, in the
    //! transverse plane, of the projectile's hard quark at proj_q
    SpatialVec sample_targ_hard_quark_position(
        std::shared_ptr<Nucleon> targ, const SpatialVec &targ_n,
        const SpatialVec &proj_q, real BG);
    void Set_hard_parton_momentum(std::vector<double> &HardMomandPosProj,
                                  std::vector<double> &HardMomandPosTarg);
    void Set_hard_collisions_Pos(std::vector<double> &HardPosProj);

    void output_spectators(std::string filename);
    void prepareParticipantList();
    std::vector<std::vector<real>> getParticipantList() {
        prepareParticipantList();
        return(participantList_);
    }
    void outputParticipants(std::string filename);
    std::vector<std::vector<real>> get_QCD_strings_output_list() {
        prepare_output_QCD_strings();
        return(QCD_string_output_arr_);
    }

    //! The string output of this event with a hard process's energy taken
    //! from the two nucleons of binary collision collision_idx (an index of
    //! get_collision_information()), for a second hydro run on the same
    //! event: a modified copy, the event itself is not changed, and every
    //! string keeps the event's centre of mass.  The event's strings must
    //! have been produced (perform_string_production) with valence quarks.
    //!   mode 1: each nucleon's E_proj / E_targ [GeV] is taken from one of
    //!           its valence quarks that made strings (E_q >= E + m_q; picked
    //!           by a hash of the candidates' rapidities, no random number is
    //!           drawn), or, if none has enough, from several, highest E_q
    //!           first, up to 90 % of each.  The quark's rapidity is lowered
    //!           to +-acosh((E_q - E)/m_q) and all its strings move with it,
    //!           re-evolved as produced, so the strings deposit exactly E
    //!           less unless a string's rapidity loss gets capped.
    //! info, HARD_SUBTRACTION_INFO entries:
    //!   0/1 E requested proj/targ, 2/3 E taken, 4/5 quarks used,
    //!   6/7 rapidity shift of the (first) quark, 8/9 its energy E_q,
    //!   10 strings changed, 11 strings whose rapidity loss got capped,
    //!   12 strings evolve_QCD_string did not reproduce before the change,
    //!   13/14 E of the output before/after, 15/16 P_z before/after,
    //!   17 strings between the hard pair
    static constexpr int HARD_SUBTRACTION_INFO = 18;
    std::vector<std::vector<real>> get_QCD_strings_output_list_hard_subtracted(
        const int collision_idx, const real E_proj, const real E_targ,
        const int mode, std::vector<real> &info);

    real get_sig_eff(const real siginNN);
};

}

#endif   // SRC_GLAUBER_H_

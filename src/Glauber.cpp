// Copyright @ Chun Shen 2018

#include "Glauber.h"
#include "data_structs.h"
#include "PhysConsts.h"
#include <string>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <memory>
#include <cstdint>
#include <cstring>
#include <numeric>
#include <stdexcept>
#include <tuple>

using std::cout;
using std::endl;
using std::shared_ptr;

namespace MCGlb {

Glauber::Glauber(const MCGlb::Parameters &param_in,
                 shared_ptr<RandomUtil::Random> ran_gen) :
    parameter_list(param_in) {
    sample_valence_quark = false;
    fluct_Nstrings_per_NN_collision_ = false;
    remnant_energy_loss_fraction_ = (
            parameter_list.get_remnant_energy_loss_fraction());
    if (parameter_list.get_use_quarks() > 0) {
        sample_valence_quark = true;
        fluct_Nstrings_per_NN_collision_ = (
                parameter_list.get_fluct_Nstrings_per_NN_collision());
        if (!parameter_list.get_cached_tabels()) {
            system_status_ = std::system(
                    "rm -fr tables/proton_valence_quark_samples*");
            system_status_ = std::system(
                    "rm -fr tables/neutron_valence_quark_samples*");
        }
    }

    real d_min = parameter_list.get_d_min();

    int N_sea_partons = parameter_list.get_N_sea_partons();

    bool deformed = true;
    bool nucleonConfFromFile = parameter_list.nucleon_configuration_from_file();
    projectile = std::unique_ptr<Nucleus>(
            new Nucleus(parameter_list.get_projectle_nucleus_name(), ran_gen,
                        sample_valence_quark, parameter_list.get_BG(),
                        d_min, deformed, nucleonConfFromFile, N_sea_partons));
    int resetProjWS = static_cast<int>(
                            parameter_list.getParam("resetProjWS", 0.0));
    if (resetProjWS != 0) {
        auto defaultWSParam = projectile->get_woods_saxon_parameters();
        real WS_rho = parameter_list.getParam(
                                    "ProjWS_rho0", defaultWSParam[0]);
        real WS_w = parameter_list.getParam("ProjWS_w", defaultWSParam[1]);
        real WS_R = parameter_list.getParam("ProjWS_R", defaultWSParam[2]);
        real WS_a = parameter_list.getParam("ProjWS_a", defaultWSParam[3]);
        real WS_beta2 = parameter_list.getParam(
                                    "ProjWS_beta2", defaultWSParam[4]);
        real WS_beta3 = parameter_list.getParam(
                                    "ProjWS_beta3", defaultWSParam[5]);
        real WS_beta4 = parameter_list.getParam(
                                    "ProjWS_beta4", defaultWSParam[6]);
        real WS_gamma = parameter_list.getParam(
                                    "ProjWS_gamma", defaultWSParam[7]);
        projectile->setWoodsSaxonParameters(
            WS_rho, WS_w, WS_R, WS_a, WS_beta2, WS_beta3, WS_beta4, WS_gamma);
    }

    target = std::unique_ptr<Nucleus>(
            new Nucleus(parameter_list.get_target_nucleus_name(), ran_gen,
                        sample_valence_quark, parameter_list.get_BG(),
                        d_min, deformed, nucleonConfFromFile, N_sea_partons));
    int resetTargWS = static_cast<int>(
                            parameter_list.getParam("resetTargWS", 0.0));
    if (resetTargWS != 0) {
        auto defaultWSParam = target->get_woods_saxon_parameters();
        real WS_rho, WS_w, WS_R, WS_a, WS_beta2, WS_beta3, WS_beta4, WS_gamma;
        if (parameter_list.get_target_nucleus_name()
                == parameter_list.get_projectle_nucleus_name()) {
            // if it is a symmetric collision, then use the projectile WS param
            WS_rho = parameter_list.getParam(
                                   "ProjWS_rho0", defaultWSParam[0]);
            WS_w = parameter_list.getParam("ProjWS_w", defaultWSParam[1]);
            WS_R = parameter_list.getParam("ProjWS_R", defaultWSParam[2]);
            WS_a = parameter_list.getParam("ProjWS_a", defaultWSParam[3]);
            WS_beta2 = parameter_list.getParam(
                                   "ProjWS_beta2", defaultWSParam[4]);
            WS_beta3 = parameter_list.getParam(
                                   "ProjWS_beta3", defaultWSParam[5]);
            WS_beta4 = parameter_list.getParam(
                                   "ProjWS_beta4", defaultWSParam[6]);
            WS_gamma = parameter_list.getParam(
                                   "ProjWS_gamma", defaultWSParam[7]);
        } else {
            // for asymmetric collisions, set with target WS params
            WS_rho = parameter_list.getParam(
                                   "TargWS_rho0", defaultWSParam[0]);
            WS_w = parameter_list.getParam("TargWS_w", defaultWSParam[1]);
            WS_R = parameter_list.getParam("TargWS_R", defaultWSParam[2]);
            WS_a = parameter_list.getParam("TargWS_a", defaultWSParam[3]);
            WS_beta2 = parameter_list.getParam(
                                   "TargWS_beta2", defaultWSParam[4]);
            WS_beta3 = parameter_list.getParam(
                                   "TargWS_beta3", defaultWSParam[5]);
            WS_beta4 = parameter_list.getParam(
                                   "TargWS_beta4", defaultWSParam[6]);
            WS_gamma = parameter_list.getParam(
                                   "TargWS_gamma", defaultWSParam[7]);
        }
        target->setWoodsSaxonParameters(
            WS_rho, WS_w, WS_R, WS_a, WS_beta2, WS_beta3, WS_beta4, WS_gamma);
    }
    if (nucleonConfFromFile) {
        projectile->setLightNucleusOption(
                        parameter_list.getLightNucleusOption());
        target->setLightNucleusOption(
                        parameter_list.getLightNucleusOption());
    }

    if (sample_valence_quark) {
        projectile->set_valence_quark_Q2(parameter_list.get_quarks_Q2());
        target->set_valence_quark_Q2(parameter_list.get_quarks_Q2());
    }
    ran_gen_ptr_ = ran_gen;

    yloss_param_slope = parameter_list.get_yloss_param_slope();
    real alpha1 = parameter_list.get_yloss_param_alpha1();
    real alpha2 = parameter_list.get_yloss_param_alpha2();
    real alpha = alpha2/alpha1;
    yloss_param_a = alpha/(1. - alpha);
    yloss_param_b = alpha2/yloss_param_a;

    ybeam = acosh(parameter_list.get_roots()/(2.*PhysConsts::MProton));

    real siginNN = compute_NN_inelastic_cross_section(
                                            parameter_list.get_roots());
    sigma_eff_ = get_sig_eff(siginNN);
    cout << "sqrt{s} = " << parameter_list.get_roots() << " GeV, "
         << "siginNN = " << siginNN << " mb" << endl;
}

void Glauber::make_nuclei() {
    projectile->generate_nucleus_3d_configuration();
    target->generate_nucleus_3d_configuration();
    projectile->accelerate_nucleus(parameter_list.get_roots(), 1);
    target->accelerate_nucleus(parameter_list.get_roots(), -1);

    // sample impact parameters
    auto b_max = parameter_list.get_b_max();
    auto b_min = parameter_list.get_b_min();
    impact_b = sqrt(b_min*b_min +
            (b_max*b_max - b_min*b_min)*ran_gen_ptr_->rand_uniform());
    SpatialVec proj_shift = {0., impact_b/2., 0.,
                             -projectile->get_z_max() - 1e-15};
    projectile->shift_nucleus(proj_shift);
    SpatialVec targ_shift = {0., -impact_b/2., 0.,
                             -target->get_z_min() + 1e-15};
    target->shift_nucleus(targ_shift);
    //projectile->output_nucleon_positions("projectile.dat");
    //target->output_nucleon_positions("target.dat");
}


real Glauber::get_nucleus_density(double t, double x,
                                  double y, double z, int direction,
                                  std::unique_ptr<Nucleus> &nucleus_ptr) {
    auto nuclues_rapidity = nucleus_ptr->get_beam_rapidity(
                                    parameter_list.get_roots(), direction);
    const real width_xy = 0.5;      // fm
    const real width_z = width_xy/cosh(nuclues_rapidity);   // fm
    const real prefactor = 1./(width_xy*width_xy*width_z)/pow(M_PI*2, 1.5);

    auto nucleon_list = nucleus_ptr->get_nucleon_list();
    real nucleon_density = 0.0;
    for (auto &inucleon: (*nucleon_list)) {
        auto nucl_x = inucleon->get_x();
        auto dis_square_xy = ((x - nucl_x[1])*(x - nucl_x[1]) +
                              (y - nucl_x[2])*(y - nucl_x[2]));
        if (dis_square_xy > width_xy*width_xy*200.) continue;

        auto z_shifted = z - tanh(nuclues_rapidity)*t;
        auto dis_square_z = (z_shifted - nucl_x[3])*(z_shifted - nucl_x[3]);
        if (dis_square_z > width_z*width_z*100.) continue;

        nucleon_density += exp(- dis_square_xy/(2.*width_xy*width_xy)
                               - dis_square_z/(2.*width_z*width_z));
    }
    return (nucleon_density*prefactor);
}


real Glauber::get_proj_nucleon_density(double t, double x,
                                       double y, double z) {
    return get_nucleus_density(t, x, y, z, 1, projectile);
}


real Glauber::get_targ_nucleon_density(double t, double x,
                                       double y, double z) {
    return get_nucleus_density(t, x, y, z, -1, target);
}


real Glauber::get_nucleon_density(double t, double x, double y, double z) {
    auto nucleon_density = (get_proj_nucleon_density(t, x, y, z)
                            + get_targ_nucleon_density(t, x, y, z));
    return (nucleon_density);
}


std::vector<double>  Glauber::get_all_proj_nucleon_z() {
    Proj_nucleonz_.clear();
    auto nucleon_list = projectile->get_nucleon_list();
    for (auto &inucleon: (*nucleon_list)) {
        auto nucl_x = inucleon->get_x();
        Proj_nucleonz_.push_back(nucl_x[3]);
    }
    return (Proj_nucleonz_);
}


std::vector<double>  Glauber::get_all_targ_nucleon_z() {
    Targ_nucleonz_.clear();
    auto nucleon_list = target->get_nucleon_list();
    for (auto &inucleon: (*nucleon_list)) {
        auto nucl_x = inucleon->get_x();
        Targ_nucleonz_.push_back(nucl_x[3]);
    }
    return (Targ_nucleonz_);
}


int Glauber::make_collision_schedule() {
    collision_schedule.clear();
    auto proj_nucleon_list = projectile->get_nucleon_list();
    auto targ_nucleon_list = target->get_nucleon_list();
    for (auto &iproj: (*proj_nucleon_list)) {
        auto proj_x = iproj->get_x();
        for (auto &itarg: (*targ_nucleon_list)) {
            auto targ_x = itarg->get_x();
            auto dij = (  (targ_x[1] - proj_x[1])*(targ_x[1] - proj_x[1])
                        + (targ_x[2] - proj_x[2])*(targ_x[2] - proj_x[2]));
            if (hit(dij)) {
                create_a_collision_event(iproj, itarg);
                projectile->add_a_participant(iproj);
                target->add_a_participant(itarg);
                iproj->set_wounded(true);
                itarg->set_wounded(true);
                iproj->add_collide_nucleon(itarg);
                itarg->add_collide_nucleon(iproj);
            }
        }
    }

    // copy the collision_schedule information to outside
    collision_schedule_list_.clear();
    for (auto &it: collision_schedule) {
        // collision list is time ordered
        collision_schedule_list_.push_back(*it);
    }
    get_collision_information();
    return(collision_schedule.size());
}

int Glauber::make_collision_schedule_second_stage() {
    collision_schedule.clear();
    auto proj_nucleon_list = projectile->get_nucleon_list();
    auto targ_nucleon_list = target->get_nucleon_list();
    for (auto &iproj: (*proj_nucleon_list)) {
        for (auto &itarg: (*targ_nucleon_list)) {
            create_a_collision_event(iproj, itarg);
            projectile->add_a_participant(iproj);
            target->add_a_participant(itarg);
            iproj->set_wounded(true);
            itarg->set_wounded(true);
            iproj->add_collide_nucleon(itarg);
            itarg->add_collide_nucleon(iproj);
        }
    }

    // copy the collision_schedule information to outside
    collision_schedule_list_.clear();
    for (auto &it: collision_schedule) {
        // collision list is time ordered
        collision_schedule_list_.push_back(*it);
    }
    get_collision_information();
    return(collision_schedule.size());
}


int Glauber::get_Npart() const {
    int Npart = (projectile->get_number_of_wounded_nucleons()
                  + target->get_number_of_wounded_nucleons());
    return(Npart);
}


bool Glauber::hit(real d2) const {
    //real G = 0.92;  // from Glassando
    //return(ran_gen_ptr_->rand_uniform() < G*exp(-G*d2/d2_in));
    const real T_nn = (exp(-d2/(4.*nucleon_width_*nucleon_width_))
                       /(4.*M_PI*nucleon_width_*nucleon_width_));
    const real hit_treshold = 1. - exp(-sigma_eff_*T_nn);
    return (ran_gen_ptr_->rand_uniform() < hit_treshold);
}


void Glauber::create_a_collision_event(shared_ptr<Nucleon> proj,
                                       shared_ptr<Nucleon> targ) {
    real t_coll, z_coll;
    auto x1 = proj->get_x();
    auto p1 = proj->get_p();
    auto x2 = targ->get_x();
    auto p2 = targ->get_p();
    auto v1 = p1[3]/p1[0];
    auto v2 = p2[3]/p2[0];
    auto collided = get_collision_point(x1[0], x1[3], v1, x2[3], v2,
                                        t_coll, z_coll);
    if (collided) {
        SpatialVec x_coll = {t_coll, (x1[1] + x2[1])/2.,
                             (x1[2] + x2[2])/2., z_coll};
        shared_ptr<CollisionEvent> event_ptr(
                                    new CollisionEvent(x_coll, proj, targ));
        if (proj->is_connected_with(targ)) {
            auto form_N_strings = proj->get_number_of_connections(targ);
            event_ptr->set_produced_n_strings(form_N_strings);
        }
        collision_schedule.insert(event_ptr);
    }
}


bool Glauber::get_collision_point(real t, real z1, real v1, real z2, real v2,
                                  real &t_coll, real &z_coll) const {
    bool collided = false;
    real delta_t = (z2 - z1)/(v1 - v2);
    if (delta_t > 0.) {
        collided = true;
        t_coll = t  + delta_t;
        z_coll = z1 + v1*delta_t;
    } else {
        t_coll = -1.;
        z_coll = -1.;
    }
    return(collided);
}


real Glauber::compute_NN_inelastic_cross_section(real ecm) const {
    real s = ecm*ecm;
    real sigma_NN_total = 44.4 - 2.9*log(s) + 0.33*log(s)*log(s);
    real sigma_NN_inel  = (sigma_NN_total
                           - (11.4 - 1.52*log(s) + 0.13*log(s)*log(s)));
    return(sigma_NN_inel);
}


void Glauber::Set_hard_parton_momentum(std::vector<double> &HardMomandPosProj,
                                       std::vector<double> &HardMomandPosTarg) {
    HardPartonPosAndMomProj_.clear();
    HardPartonPosAndMomTarg_.clear();
    for (unsigned int i=0; i<HardMomandPosProj.size(); i++) {
        HardPartonPosAndMomProj_.push_back(HardMomandPosProj[i]);
        HardPartonPosAndMomTarg_.push_back(HardMomandPosTarg[i]);
    }
}


real Glauber::GG_probality(real x, real Lambda, real Omega, real Sigma0) const {
    return( 1./Lambda * x/Lambda / (x/Lambda + Sigma0) *
            exp(-1 * pow(x/Lambda/Sigma0 - 1., 2.) / Omega/Omega));
}

void Glauber::Set_hard_collisions_Pos(std::vector<double> &HardPosProj) {
    HardPartonPos_.clear();
    for (unsigned int i=0; i<HardPosProj.size(); i++) {
        HardPartonPos_.push_back(HardPosProj[i]);
    }
}

std::vector<double>  Glauber::OutputquarkPosProj() {
    auto binary_collision_x = HardPartonPos_[1];
    auto binary_collision_y = HardPartonPos_[2];
    std::vector<double> quark_xvec;
    // pick up the colliding nucleon pair generated hard partons
    for (auto &it: collision_schedule_list_) {
        // collision list is time ordered
        auto xvec = it.get_collision_position();
        if ((abs(xvec[1] - binary_collision_x) < 1.e-5) && 
            (abs(xvec[2] - binary_collision_y) < 1.e-5)) {
             // Pick up one valence quark,
             auto nucleon_collided = it.get_proj_nucleon_ptr().lock();
             if (sample_valence_quark) {
                 quark_xvec = nucleon_collided->output_quark_pos();
                 break;
            }
        }
    }
    return(quark_xvec);
}

std::vector<double>  Glauber::OutputquarkPosTarg() {
    auto binary_collision_x = HardPartonPos_[1];
    auto binary_collision_y = HardPartonPos_[2];
    std::vector<double> quark_xvec;
    // pick up the colliding nucleon pair generated hard partons
    for (auto &it: collision_schedule_list_) {
        // collision list is time ordered
        auto xvec = it.get_collision_position();
        if ((abs(xvec[1] - binary_collision_x) < 1.e-5) && 
            (abs(xvec[2] - binary_collision_y) < 1.e-5)) {
             // Pick up one valence quark,
             auto nucleon_collided = it.get_targ_nucleon_ptr().lock();
             if (sample_valence_quark) {
                 quark_xvec = nucleon_collided->output_quark_pos();
                 break;
            }
        }
    }
    return(quark_xvec);
}

void Glauber::Pick_and_subtract_hard_parton_momentum() {
    real ecm = parameter_list.get_roots();
    const real mN = PhysConsts::MProton;
    real nucleonPz = mN*sinh(ybeam);
    // Positions and Momentum for the leading hard partons.
    auto binary_collision_x = HardPartonPosAndMomProj_[1];
    auto binary_collision_y = HardPartonPosAndMomProj_[2];
    // pick up the colliding nucleon pair generated hard partons
    for (auto &it: collision_schedule_list_) {
        // collision list is time ordered
        auto xvec = it.get_collision_position();
        if ((std::abs(xvec[1] - binary_collision_x) < 1.e-5) &&
            (std::abs(xvec[2] - binary_collision_y) < 1.e-5)) {
            // Pick up one valence quark,
            // substract the four momentum from this valence quark
            auto proj_collided = it.get_proj_nucleon_ptr().lock();
            auto targ_collided = it.get_targ_nucleon_ptr().lock();
            proj_collided->set_hard_collided(true);
            targ_collided->set_hard_collided(true);
            if (!parameter_list.subtract_hard_momentum()) {
                HardPartonPosAndMomProj_[4] = 0.0;
                HardPartonPosAndMomProj_[5] = 0.0;
                HardPartonPosAndMomProj_[6] = 0.0;
                HardPartonPosAndMomProj_[7] = 0.0;

                HardPartonPosAndMomTarg_[4] = 0.0;
                HardPartonPosAndMomTarg_[5] = 0.0;
                HardPartonPosAndMomTarg_[6] = 0.0;
                HardPartonPosAndMomTarg_[7] = 0.0;
            }
            MomentumVec HardPartonMomProj_ = { HardPartonPosAndMomProj_[4], 
                                               HardPartonPosAndMomProj_[5],
                                               HardPartonPosAndMomProj_[6],
                                               HardPartonPosAndMomProj_[7] };
            MomentumVec HardPartonMomTarg_ = { HardPartonPosAndMomTarg_[4], 
                                               HardPartonPosAndMomTarg_[5],
                                               HardPartonPosAndMomTarg_[6],
                                               HardPartonPosAndMomTarg_[7] };
            // Pick up the valence quark
            if (sample_valence_quark) {
                if (HardPartonMomProj_[0] > ecm/2.1
                       || HardPartonMomTarg_[0] > ecm/2.1) {
                    std::cout << "Sorry for the too large hard collision energy,";
                    std::cout << " Bye ~~" << std::endl;
                    exit(-1);
                }
                std::shared_ptr<Quark> proj_q;
                std::shared_ptr<Quark> targ_q;
                MomentumVec p_q;
                int do_resample_proj = 1;
                SpatialVec pos_q_Proj, pos_q_Targ;
                while (do_resample_proj == 1) {
                    proj_q = proj_collided->get_a_valence_quark_sub_mom(
                                                   HardPartonMomProj_[0]);
                    p_q = proj_q->get_p();
                    if (p_q[0] <= HardPartonMomProj_[0]) {
                        // resample the valence quark requiring
                        // one parton with x >= x_hard
                        real x_hard = std::abs(HardPartonMomProj_[3]/nucleonPz);
                        std::vector<double> xvec_q = (
                                       proj_collided->output_quark_pos());

                        std::vector<real> xQuark;
                        projectile->resample_quark_momentum_fraction(
                                xQuark,
                                proj_collided->get_electric_charge(),
                                ecm, x_hard);
                        proj_collided->resample_valence_quarks(ecm, 1, xvec_q,
                                                               xQuark);
                        proj_collided->readd_soft_parton_ball(
                               ecm, 1, xvec_q, parameter_list.get_BG(),
                               proj_collided->get_p(),
                               proj_collided->get_quark_list());
                        //std::cout << " re-sample the valence quark in proj." <<std::endl;
                    } else {
                        proj_q->set_subtracted(true);
                        auto proj_q_pos = proj_q->get_x();
                        auto proj_n_pos = proj_collided->get_x();
                         pos_q_Proj = {0.0,
                                       proj_q_pos[1] + proj_n_pos[1],
                                       proj_q_pos[2] + proj_n_pos[2],
                                       proj_q_pos[3] + proj_n_pos[3]};
                         set_Proj_hot_spot_x(proj_q->get_x());
                         do_resample_proj = 0;
                     }
                }
                proj_collided->erase_one_quark();
                int do_resample_targ = 1;
                while (do_resample_targ == 1) {
                    targ_q = targ_collided->get_a_valence_quark_sub_mom(
                                                        HardPartonMomTarg_[0]);
                    p_q = targ_q->get_p();
                    if (p_q[0] <= HardPartonMomTarg_[0]) { // 0.346**2~ 0.12
                        // resample the valence quark requiring
                        // one parton with x >= x_hard
                        real x_hard = std::abs(HardPartonMomTarg_[3]/nucleonPz);
                        std::vector<double> xvec_q = targ_collided->output_quark_pos();
                        std::vector<real> xQuark;
                        target->resample_quark_momentum_fraction(
                                xQuark,
                                targ_collided->get_electric_charge(),
                                ecm, x_hard);
                        targ_collided->resample_valence_quarks(ecm, -1, xvec_q,
                                                               xQuark);
                        targ_collided->readd_soft_parton_ball(
                               ecm, -1, xvec_q, parameter_list.get_BG(),
                               targ_collided->get_p(),
                               targ_collided->get_quark_list());
                        //std::cout << " re-sample the valence quark in targ." <<std::endl;
                    } else {
                        targ_q->set_subtracted(true);
                        SpatialVec tarj_q_pos = targ_q->get_x();
                        auto tarj_n = targ_collided->get_x();
                        pos_q_Targ = {0.0,
                                      tarj_q_pos[1] + tarj_n[1],
                                      tarj_q_pos[2] + tarj_n[2],
                                      tarj_q_pos[3] + tarj_n[3]};
                        auto dis_square_q_xy = (pos_q_Targ[1] - pos_q_Proj[1])*(pos_q_Targ[1] - pos_q_Proj[1]) +
                                               (pos_q_Targ[2] - pos_q_Proj[2])*(pos_q_Targ[2] - pos_q_Proj[2]);
                        if (dis_square_q_xy > 0.25) {
                            tarj_q_pos = sample_targ_hard_quark_position(
                                targ_collided, tarj_n, pos_q_Proj,
                                parameter_list.get_BG());
                        }
                        targ_q->set_x(tarj_q_pos);
                        set_Targ_hot_spot_x(tarj_q_pos);
                        do_resample_targ = 0;
                    }
                }
                targ_collided->erase_one_quark();
                break;
            } else {
                proj_collided->substract_momentum_from_remnant(HardPartonMomProj_);
                targ_collided->substract_momentum_from_remnant(HardPartonMomTarg_);
                break;
            }
        }
    }
}


SpatialVec Glauber::sample_targ_hard_quark_position(
        std::shared_ptr<Nucleon> targ, const SpatialVec &targ_n,
        const SpatialVec &proj_q, real BG) {
    // Whole positions from the valence-quark Gaussian until one falls within
    // 0.5 fm of the projectile's quark.  Its chance per draw falls like
    // exp(-d^2/(2 sigma^2)) with the distance d between the projectile's
    // quark and the target nucleon (sigma = sqrt(BG) hbarc ~ 0.4 fm), so
    // beyond d ~ 2.5 fm this would practically never end (seen in ~1 of
    // 500-1300 O+O events with ISR).  1e8 draws take ~5 s; an event that
    // ends within them keeps its random numbers.
    const long max_draws = 100000000;
    for (long n = 0; n < max_draws; n++) {
        SpatialVec x_q = targ->resample_valence_quark_position(BG);
        const real dx = x_q[1] + targ_n[1] - proj_q[1];
        const real dy = x_q[2] + targ_n[2] - proj_q[2];
        if (dx*dx + dy*dy <= 0.25) return(x_q);
    }
    // The same distribution, sampled directly: uniform in the 0.5 fm disk
    // around the projectile's quark, kept with the Gaussian's weight relative
    // to the disk's point nearest the nucleon (<= 1).  z as in the draws.
    const real sigma = sqrt(BG)*PhysConsts::HBARC;
    const real R = 0.5;
    const real cx = proj_q[1] - targ_n[1];   // disk centre, from the nucleon
    const real cy = proj_q[2] - targ_n[2];
    const real d_min = std::max(0., sqrt(cx*cx + cy*cy) - R);
    real x, y;
    do {
        const real r = R*sqrt(ran_gen_ptr_->rand_uniform());
        const real phi = 2.*M_PI*ran_gen_ptr_->rand_uniform();
        x = cx + r*cos(phi);
        y = cy + r*sin(phi);
    } while (ran_gen_ptr_->rand_uniform()
             >= exp(-(x*x + y*y - d_min*d_min)/(2.*sigma*sigma)));
    const real z = ran_gen_ptr_->rand_normal(0., sigma);
    std::cout << "[3dMCGlauber] the target's hard quark was placed by direct "
              << "sampling (" << max_draws << " draws missed; projectile "
              << "quark " << sqrt(cx*cx + cy*cy) << " fm from the target "
              << "nucleon)" << std::endl;
    return(SpatialVec{0.0, x, y, z});
}


int Glauber::decide_produce_string_num(
                                shared_ptr<CollisionEvent> event_ptr) const {
    int form_n_string = 0;
    auto proj = event_ptr->get_proj_nucleon_ptr().lock();
    auto targ = event_ptr->get_targ_nucleon_ptr().lock();
    int minimum_allowed_connections = 1;
    if (sample_valence_quark && fluct_Nstrings_per_NN_collision_) {
        // assume P(N) = 1/e^N distribution
        int N = std::min(proj->get_number_of_quarks(),
                         targ->get_number_of_quarks());
        do {
            auto rand = ran_gen_ptr_->rand_uniform();
            minimum_allowed_connections = static_cast<int>(-log(1 - rand)) + 1;
        } while (minimum_allowed_connections > N);
    }

    if (   proj->get_number_of_connections() < minimum_allowed_connections
        && targ->get_number_of_connections() < minimum_allowed_connections) {
        form_n_string = std::min(
            minimum_allowed_connections - proj->get_number_of_connections(),
            minimum_allowed_connections - targ->get_number_of_connections());
    } else if (
           proj->get_number_of_connections() < minimum_allowed_connections
        || targ->get_number_of_connections() < minimum_allowed_connections) {
        form_n_string = 1;
    } else {
        int n_connects = (  proj->get_number_of_connections()
                          + targ->get_number_of_connections());
        real shadowing_factor = parameter_list.get_shadowing_factor();
        real production_prob = (
            (1. - shadowing_factor)
            *exp(-shadowing_factor*(n_connects
                                    - 2.*minimum_allowed_connections)));
        if (ran_gen_ptr_->rand_uniform() < production_prob) {
            form_n_string = 1;
        }
        if (proj->is_hard_collided() && targ->is_hard_collided()) {
            form_n_string = 0;
        }
    }
    return(form_n_string);
}


void Glauber::sample_valence_quarks() {
    if (sample_valence_quark) {
        projectile->sample_valence_quarks_inside_nucleons(
                                    parameter_list.get_roots(), 1);
        target->sample_valence_quarks_inside_nucleons(
                                    parameter_list.get_roots(), -1);
        projectile->add_soft_parton_ball(parameter_list.get_roots(), 1);
        target->add_soft_parton_ball(parameter_list.get_roots(), -1);
    }
}


int Glauber::decide_QCD_strings_production() {
    sample_valence_quarks();
    std::vector<shared_ptr<CollisionEvent>> collision_list;
    for (auto &it: collision_schedule)
        collision_list.push_back(it);  // collision list is time ordered

    const auto QCD_string_production_mode =
                            parameter_list.get_QCD_string_production_mode();
    if (QCD_string_production_mode == 1) {
        // randomly ordered strings
        std::shuffle(collision_list.begin(), collision_list.end(),
                     *ran_gen_ptr_->getRanGenerator());
    } else if (QCD_string_production_mode == 2) {
        // anti-time ordered strings
        std::reverse(collision_list.begin(), collision_list.end());
    }

    int number_of_strings = 0;
    bool finished = false;
    while (!finished) {
        finished = true;
        for (auto &ievent: collision_list) {
            auto form_N_strings = decide_produce_string_num(ievent);
            auto proj = ievent->get_proj_nucleon_ptr();
            auto targ = ievent->get_targ_nucleon_ptr();
            if (form_N_strings > 0) {
                number_of_strings += form_N_strings;
                proj.lock()->add_connected_nucleon(targ);
                proj.lock()->add_num_connections(form_N_strings);
                targ.lock()->add_connected_nucleon(proj);
                targ.lock()->add_num_connections(form_N_strings);
                ievent->set_produced_n_strings(form_N_strings);
            } else {
                if (   proj.lock()->get_number_of_connections() == 0
                    || targ.lock()->get_number_of_connections() == 0)
                        finished = false;
            }
        }
    }
    return(number_of_strings);
}


int Glauber::decide_QCD_strings_production_second_stage() {
    std::vector<shared_ptr<CollisionEvent>> collision_list;
    collision_list.clear();
    for (auto &it: collision_schedule)
        collision_list.push_back(it);  // collision list is time ordered

    const auto QCD_string_production_mode =
                            parameter_list.get_QCD_string_production_mode();
    if (QCD_string_production_mode == 1) {
        // randomly ordered strings
        std::shuffle(collision_list.begin(), collision_list.end(),
                     *ran_gen_ptr_->getRanGenerator());
    } else if (QCD_string_production_mode == 2) {
        // anti-time ordered strings
        std::reverse(collision_list.begin(), collision_list.end());
    }

    int number_of_strings = 0;
    bool finished = false;
    while (!finished) {
        finished = true;
        for (auto &ievent: collision_list) {
            auto form_N_strings = decide_produce_string_num(ievent);
            auto proj = ievent->get_proj_nucleon_ptr();
            auto targ = ievent->get_targ_nucleon_ptr();
            ievent->set_produced_n_strings(form_N_strings);
            if (form_N_strings > 0) {
                number_of_strings += form_N_strings;
                proj.lock()->add_connected_nucleon(targ);
                proj.lock()->add_num_connections(form_N_strings);
                targ.lock()->add_connected_nucleon(proj);
                targ.lock()->add_num_connections(form_N_strings);
                //ievent->set_produced_n_strings(form_N_strings);
            } else {
                if (   proj.lock()->get_number_of_connections() == 0
                    || targ.lock()->get_number_of_connections() == 0)
                        finished = false;
            }
        }
    }
    return(number_of_strings);
}

//! This function propagate individual nucleon inside the nucleus by dt
void Glauber::propagate_nuclei(real dt) {
    for (auto &n_i: (*projectile->get_nucleon_list()))
        propagate_nucleon(n_i, dt);
    for (auto &n_i: (*target->get_nucleon_list()))
        propagate_nucleon(n_i, dt);
}


void Glauber::propagate_nucleon(shared_ptr<Nucleon> n_i, real dt) {
    auto x_vec = n_i->get_x();
    auto p_vec = n_i->get_p();
    for (int i = 0; i < 4; i++) {
        real v_i = p_vec[i]/p_vec[0];
        x_vec[i] += dt*v_i;
    }
    n_i->set_x(x_vec);
}


void Glauber::update_momentum_quark(shared_ptr<Quark> q_i, real y_shift) {
    auto pvec = q_i->get_p();
    auto mass = q_i->get_mass();
    auto y_i = q_i->get_rapidity();
    auto y_f = y_i + y_shift;
    q_i->set_rapidity(y_f);
    pvec[0] = mass*cosh(y_f);
    pvec[3] = mass*sinh(y_f);
    q_i->set_p(pvec);
}


void Glauber::update_momentum(shared_ptr<Nucleon> n_i, real y_shift) {
    auto pvec = n_i->get_p();
    auto y_i = n_i->get_rapidity();
    auto y_f = y_i + y_shift;
    pvec[0] = PhysConsts::MProton*cosh(y_f);
    pvec[3] = PhysConsts::MProton*sinh(y_f);
    n_i->set_p(pvec);
}


real Glauber::get_tau_form(const int string_evolution_mode) const {
    auto tau_form_mean = parameter_list.get_tau_form_mean();
    real tau_form = tau_form_mean;      // [fm]
    if (string_evolution_mode == 2) {
        // tau_form fluctuates with a gamma distribution
        tau_form = tau_form_mean/ran_gen_ptr_->rand_gamma_dis();
        tau_form = std::min(std::max(0.1, tau_form), 6*tau_form_mean);
    }
    return(tau_form);
}


void Glauber::get_tau_form_and_moversigma(const int string_evolution_mode,
                                          const real y_in_lrf,
                                          real &tau_form, real &m_over_sigma,
                                          real &y_loss) {
    tau_form = get_tau_form(string_evolution_mode);      // [fm]
    m_over_sigma = 1.0;  // [fm]
    y_loss = 0.;
    real eps = 1e-16;
    if (string_evolution_mode == 1) {
        // fixed rapidity loss
        y_loss = (
            acosh(tau_form*tau_form/(2.*m_over_sigma*m_over_sigma) + 1.));
        // maximum 90% of y_in_lrf
        y_loss = std::min(y_loss, y_in_lrf*0.9);
        m_over_sigma = tau_form/std::max(eps, sqrt(2.*(cosh(y_loss) - 1.)));
    } else if (string_evolution_mode == 2) {
        // both tau_form and sigma fluctuate
        y_loss = sample_rapidity_loss_shell(y_in_lrf);
        m_over_sigma = tau_form/std::max(eps, sqrt(2.*(cosh(y_loss) - 1.)));
    } else if (string_evolution_mode == 3) {
        // only tau_form fluctuates
        y_loss = sample_rapidity_loss_shell(y_in_lrf);
        tau_form = m_over_sigma*sqrt(2.*(cosh(y_loss) - 1.));
    } else if (string_evolution_mode == 4) {
        // only m_over_sigma fluctuates
        y_loss = sample_rapidity_loss_shell(y_in_lrf);
        m_over_sigma = tau_form/std::max(eps, sqrt(2.*(cosh(y_loss) - 1.)));
    } else if (string_evolution_mode == -4) {
        // only m_over_sigma fluctuates for Beam Remnants
        y_loss = (sample_rapidity_loss_shell(y_in_lrf)
                  *remnant_energy_loss_fraction_);
        m_over_sigma = tau_form/std::max(eps, sqrt(2.*(cosh(y_loss) - 1.)));
    }
}

// add flag if baryon used to proj and target
// - check it - only put in baryon if not used yet, then set it to used
int Glauber::perform_string_production() {
    QCD_string_list.clear();
    remnant_string_list_.clear();
    const auto string_evolution_mode = (
                    parameter_list.get_QCD_string_evolution_mode());
    const auto baryon_junctions = parameter_list.get_baryon_junctions();

    // sqrt(parameter_list.get_roots()); // ~s^{-1/4} 
    real lambdaB = parameter_list.get_lambdaB();
    lambdaB = std::min(1., lambdaB);
    real lambdaBs = parameter_list.get_lambdaBs();
    lambdaBs = std::min(1., lambdaBs);

    real t_current = 0.0;
    int number_of_collided_events = 0;
    while (collision_schedule.size() > 0) {
        auto first_event = (*collision_schedule.begin());
        if (!first_event->is_valid()) {
            collision_schedule.erase((*collision_schedule.begin()));
            continue;
        }
        number_of_collided_events++;
        real dt = first_event->get_collision_time() - t_current;
        assert(dt > 0.);
        t_current = first_event->get_collision_time();
        propagate_nuclei(dt);

        auto proj = first_event->get_proj_nucleon_ptr().lock();
        auto targ = first_event->get_targ_nucleon_ptr().lock();
        if (proj->is_hard_collided() && !proj->nucleon_is_subtracted()) {
            MomentumVec HardPartonMomProj_ = { HardPartonPosAndMomProj_[4],
                                               HardPartonPosAndMomProj_[5],
                                               HardPartonPosAndMomProj_[6],
                                               HardPartonPosAndMomProj_[7] };
            proj->substract_momentum_from_remnant(HardPartonMomProj_);
            proj->set_hard_subtracted(true);
            //std::cout << "Subtract four momentum from picked up nucleon in proj." << std::endl;
        }
        if (targ->is_hard_collided() && !targ->nucleon_is_subtracted()) {
            MomentumVec HardPartonMomTarg_ = { HardPartonPosAndMomTarg_[4],
                                               HardPartonPosAndMomTarg_[5],
                                               HardPartonPosAndMomTarg_[6],
                                               HardPartonPosAndMomTarg_[7] };
            targ->substract_momentum_from_remnant(HardPartonMomTarg_);
            targ->set_hard_subtracted(true);
            //std::cout << "Subtract four momentum from picked up nucleon in targ." << std::endl;
        }

        if (first_event->get_num_strings() == 0) {
            collision_schedule.erase((*collision_schedule.begin()));
            continue;
        }
        for (int istring = 0; istring < first_event->get_num_strings();
                istring++) {
            auto x_coll = first_event->get_collision_position();
            auto proj = first_event->get_proj_nucleon_ptr().lock();
            auto targ = first_event->get_targ_nucleon_ptr().lock();
            real y_in_lrf = std::abs(proj->get_rapidity()
                                     - targ->get_rapidity())/2.;
            std::shared_ptr<Quark> proj_q;
            std::shared_ptr<Quark> targ_q;
            if (sample_valence_quark) {
                proj_q = proj->get_a_valence_quark();
                if (proj_q->get_number_of_connections() == 1) {
                    // first time pick-up the valence quark
                    // we need to substract the valence quark energy-momentum
                    // from the nucleon remnant energy-momentum vector
                    auto p_q = proj_q->get_p();
                    proj->substract_momentum_from_remnant(p_q);
                    proj_q->set_connected(true);
                }

                targ_q = targ->get_a_valence_quark();
                if (targ_q->get_number_of_connections() == 1) {
                    // first time pick-up the valence quark
                    // we need to substract the valence quark energy-momentum
                    // from the nucleon remnant energy-momentum vector
                    auto p_q = targ_q->get_p();
                    targ->substract_momentum_from_remnant(p_q);
                    targ_q->set_connected(true);
                }
                y_in_lrf = std::abs(  proj_q->get_rapidity()
                                    - targ_q->get_rapidity())/2.;
            }
            real tau_form = 0.5;      // [fm]
            real m_over_sigma = 1.0;  // [fm]
            real y_loss = 0.;
            get_tau_form_and_moversigma(string_evolution_mode, y_in_lrf,
                                        tau_form, m_over_sigma, y_loss);
            // set variables in case of no baryon junction transport
            bool has_baryon_left = false;
            bool has_baryon_right = false;
            if (!sample_valence_quark) {
                QCDString qcd_string(x_coll, tau_form, proj, targ, m_over_sigma,
                                     has_baryon_right, has_baryon_left);
                QCD_string_list.push_back(qcd_string);
            } else {
                auto proj_q_xvec = proj_q->get_x();
                auto targ_q_xvec = targ_q->get_x();
                SpatialVec x_coll_q = {
                    x_coll[0],
                    x_coll[1] + (proj_q_xvec[1] + targ_q_xvec[1])/2.,
                    x_coll[2] + (proj_q_xvec[2] + targ_q_xvec[2])/2.,
                    x_coll[3]};
                QCDString qcd_string(x_coll_q, tau_form, proj, targ,
                                     proj_q, targ_q, m_over_sigma,
                                     has_baryon_right, has_baryon_left);
                QCD_string_list.push_back(qcd_string);
            }
            real y_shift = y_loss;
            if (!sample_valence_quark) {
                update_momentum(proj, -y_shift);
                update_momentum(targ,  y_shift);
            } else {
                // shift the nucleon rapidity to avoid identical collision when
                // update the collision schedule
                update_momentum(proj, -1e-3*proj->get_rapidity());
                update_momentum(targ,  1e-3*targ->get_rapidity());
                update_momentum_quark(proj_q, -y_shift);
                update_momentum_quark(targ_q,  y_shift);
            }
        }
        update_collision_schedule(first_event);
        collision_schedule.erase((*collision_schedule.begin()));
    }

    // randomize the QCD_string_list and assign the baryon charge to
    // the strings
    std::vector<unsigned int> random_idx;
    const real baryonInStringProb = parameter_list.get_baryon_in_string_prob();
    unsigned int Nstrings = QCD_string_list.size();
    unsigned int Npart_proj = projectile->get_number_of_wounded_nucleons();
    unsigned int Npart_targ = target->get_number_of_wounded_nucleons();
    unsigned int total_length = Nstrings + Npart_proj + Npart_targ;
    for (unsigned int idx = 0; idx < total_length; idx++)
        random_idx.push_back(idx);
    std::shuffle(random_idx.begin(), random_idx.end(),
                 *ran_gen_ptr_->getRanGenerator());
    for (auto &idx: random_idx) {
        if (idx < Nstrings) {
            // put baryon of the projectile in the selected string
            auto proj = QCD_string_list[idx].get_proj();
            if (!proj->baryon_was_used()) {
                if (ran_gen_ptr_->rand_uniform() < baryonInStringProb) {
                    proj->set_baryon_used(true);
                    QCD_string_list[idx].set_has_baryon_right(true);
                }
            }
        } else if (idx < Nstrings + Npart_proj) {
            // put baryon of the projectile in the projectile remnant
            auto proj = projectile->get_participant(idx - Nstrings);
            auto p_i = proj->get_remnant_p();
            if (p_i[0] <= 0) continue;
            //auto mass = 0.;
            //if (std::abs(p_i[3]) < p_i[0]) {
            //    // a time-like beam remnant
            //    mass = sqrt(p_i[0]*p_i[0] - p_i[3]*p_i[3]);
            //}
            //if (!proj->baryon_was_used() && mass > 0.1) {}
            if (!proj->baryon_was_used()) {
                proj->set_baryon_used(true);
                proj->set_remnant_carry_baryon_number(true);
            }
        }
    }
    std::shuffle(random_idx.begin(), random_idx.end(),
                 *ran_gen_ptr_->getRanGenerator());
    for (auto &idx: random_idx) {
        if (idx < Nstrings) {
            // put baryon of the target in the selected string
            auto targ = QCD_string_list[idx].get_targ();
            if (!targ->baryon_was_used()) {
                if (ran_gen_ptr_->rand_uniform() < baryonInStringProb) {
                    targ->set_baryon_used(true);
                    QCD_string_list[idx].set_has_baryon_left(true);
                }
            }
        } else if (idx > Nstrings + Npart_proj - 1) {
            // put baryon of the target in the target remnant
            auto targ = target->get_participant(idx - Nstrings - Npart_proj);
            auto p_i = targ->get_remnant_p();
            if (p_i[0] <= 0) continue;
            //auto mass = 0.;
            //if (std::abs(p_i[3]) < p_i[0]) {
            //    // a time-like beam remnant
            //    mass = sqrt(p_i[0]*p_i[0] - p_i[3]*p_i[3]);
            //}
            //if (!targ->baryon_was_used() && mass > 0.1) {}
            if (!targ->baryon_was_used()) {
                targ->set_baryon_used(true);
                targ->set_remnant_carry_baryon_number(true);
            }
        }
    }

    // set baryons' rapidities
    for (auto &it: QCD_string_list) {
        it.evolve_QCD_string();
        if (!baryon_junctions) {
            // set baryon rapidities to string endpoint rapidities
            // if no junction transport is used
            it.set_final_baryon_rapidities(it.get_y_f_left(),
                                           it.get_y_f_right());
        } else {
            // sample HERE if baryon should be moved
            real y_baryon_right = 0.;
            if (it.get_has_baryon_right()) {
                if (ran_gen_ptr_->rand_uniform() < lambdaB) {
                    if (ran_gen_ptr_->rand_uniform() < lambdaBs) {
                        // y_baryon_right = sample_junction_rapidity_right(
                        //          it->get_y_i_left(), it->get_y_i_right());
                        y_baryon_right = sample_junction_rapidity_right(
                                    it.get_y_f_left(), it.get_y_f_right());
                    } else {
                        y_baryon_right = sample_junction_rapidity_uniformed(
                                    it.get_y_f_left(), it.get_y_f_right());
                    }
                } else {
                    // One should use the very initial rapidities of the
                    // colliding nucleons to determine the junction rapidity
                    // this may, however, lead to junctions lying outside the
                    // rapidity range of the final string which may cause
                    // problems in hydro and could be seen as unphysical...
                    // Another thing to discuss is the energy dependence.
                    // The cross section for baryon junctions stopping has
                    // a different root-s dependence than the usual inelastic
                    // cross section
                    // - so in principle it needs to be treated separately
                    // altogether... not sure how yet
                    y_baryon_right = it.get_y_f_right();
                }
            }
            real y_baryon_left = 0.;
            if (it.get_has_baryon_left()) {
                if (ran_gen_ptr_->rand_uniform() < lambdaB) {
                    if (ran_gen_ptr_->rand_uniform() < lambdaBs) {
                        //y_baryon_left = sample_junction_rapidity_left(
                        //          it->get_y_i_left(), it->get_y_i_right());
                        y_baryon_left = sample_junction_rapidity_left(
                                    it.get_y_f_left(), it.get_y_f_right());
                    } else {
                        y_baryon_left = sample_junction_rapidity_uniformed(
                                    it.get_y_f_left(), it.get_y_f_right());
                    }
                } else {
                    y_baryon_left = it.get_y_f_left();
                }
            }
            it.set_final_baryon_rapidities(y_baryon_left, y_baryon_right);
        }
    }


    // collision remnant is assigned to the last string which connects the
    // colliding nucleons or quarks
    for (std::vector<QCDString>::reverse_iterator it = QCD_string_list.rbegin();
            it != QCD_string_list.rend(); ++it) {
        if (sample_valence_quark) {
            // record the freeze-out space-time position for the remnant of
            // the collding nucleons at their last produced strings
            auto proj_n = it->get_proj();
            if (!proj_n->is_remnant_set()) {
                proj_n->set_remnant(true);
                auto x_frez = proj_n->get_x();
                proj_n->set_remnant_x_frez(x_frez);
            }
            auto targ_n = it->get_targ();
            if (!targ_n->is_remnant_set()) {
                targ_n->set_remnant(true);
                auto x_frez = targ_n->get_x();
                targ_n->set_remnant_x_frez(x_frez);
            }

            // set flags for quark remnants at their last connected strings
            auto proj_q = it->get_proj_q();
            if (!proj_q->is_remnant_set()) {
                proj_q->set_remnant(true);
                it->set_has_remnant_right(true);
            }
            auto targ_q = it->get_targ_q();
            if (!targ_q->is_remnant_set()) {
                targ_q->set_remnant(true);
                it->set_has_remnant_left(true);
            }
        } else {
            auto proj_n = it->get_proj();
            if (!proj_n->is_remnant_set()) {
                proj_n->set_remnant(true);
                it->set_has_remnant_right(true);
            }
            auto targ_n = it->get_targ();
            if (!targ_n->is_remnant_set()) {
                targ_n->set_remnant(true);
                it->set_has_remnant_left(true);
            }
        }
    }

    produce_remnant_strings();
    return(number_of_collided_events);
}



void Glauber::produce_remnant_strings() {
    // create strings for the beam remnants
    // clear any old remnants in case of multiple calls
    Mom_remnant_proj_.clear();
    Mom_remnant_targ_.clear();
    const auto string_evolution_mode = -4;
    real tau_form = 0.5;
    real m_over_sigma = 1.0;  // [fm]
    real y_loss = 0.;
    real mass_min = 1.e-3;
    auto proj_nucleon_list = projectile->get_nucleon_list();
    for (auto &iproj: (*proj_nucleon_list)) {
        if (iproj->is_wounded()) {
            auto x_i = iproj->get_remnant_x_frez();
            auto p_i = iproj->get_remnant_p();
            if (iproj->is_hard_collided() && iproj->nucleon_is_subtracted()) {
                for (unsigned int ip = 0; ip < p_i.size(); ip++) {
                    Mom_remnant_proj_.push_back(p_i[ip]);
                    if (parameter_list.subtract_hard_momentum()) {
                        p_i[ip] = 0.0;
                    }
                }
            }
            if (p_i[0] <= mass_min) {
                //std::cout << " No remnant." << std::endl;
                continue;
            }
            auto y_rem = ybeam;
            auto pT_i = std::sqrt(p_i[1]*p_i[1] + p_i[2]*p_i[2]);
            auto virtuality = p_i[0]*p_i[0] - pT_i*pT_i - p_i[3]*p_i[3];
            if (virtuality >= mass_min*mass_min) {
                y_rem = 0.5*log((p_i[0] + p_i[3])/(p_i[0] - p_i[3]));
            } else {
                //std::cout<< " Warning: The beam remnant of projectile is space-like vector." << std::endl;
                auto pT_temp = std::min(pT_i, std::sqrt(p_i[0]*p_i[0]-mass_min*mass_min));
                if (pT_temp < pT_i) {
                    //std::cout<< " Warning: E < pT of remnant of projectile, change pT from " << pT_i;
                    //std::cout<< " to " << pT_temp << std::endl;
                    p_i[1] = p_i[1]*pT_temp/pT_i;
                    p_i[2] = p_i[2]*pT_temp/pT_i;
                }
                //std::cout<< " Now change pz from " <<p_i[3];
                pT_i   = pT_temp;
                auto ymax = acosh(p_i[0] / std::sqrt(pT_i*pT_i + mass_min*mass_min));
                y_rem = std::min(ymax, ybeam);
                auto mT_square = p_i[0]*p_i[0]/cosh(y_rem)/cosh(y_rem);
                p_i[3] = std::sqrt(mT_square) * sinh(y_rem);
                //std::cout<< " to " << p_i[3] << std::endl;
            }
            MomentumVec targ_p_vec = {p_i[0], p_i[1], p_i[2], -p_i[3]};

            get_tau_form_and_moversigma(string_evolution_mode, y_rem,
                                        tau_form, m_over_sigma, y_loss);
            bool has_baryon_left = false;
            bool has_baryon_right = iproj->is_remnant_carry_baryon_number();
            QCDString qcd_string(x_i, tau_form, iproj, iproj,
                                 p_i, targ_p_vec, m_over_sigma,
                                 has_baryon_right, has_baryon_left);
            qcd_string.set_has_remnant_right(true);
            qcd_string.evolve_QCD_string();
            qcd_string.set_final_baryon_rapidities(0., y_rem - y_loss);
            remnant_string_list_.push_back(qcd_string);
        }
    }
    auto targ_nucleon_list = target->get_nucleon_list();
    for (auto &itarg: (*targ_nucleon_list)) {
        if (itarg->is_wounded()) {
            auto x_i = itarg->get_remnant_x_frez();
            auto p_i = itarg->get_remnant_p();
            if (itarg->is_hard_collided() && itarg->nucleon_is_subtracted()) {
                for (unsigned int ip = 0; ip < p_i.size(); ip++) {
                    Mom_remnant_targ_.push_back(p_i[ip]);
                    if (parameter_list.subtract_hard_momentum()) {
                        p_i[ip] = 0.0;
                    }
                }
            }
            if (p_i[0] <= mass_min) {
                //std::cout << " No remnant." << std::endl;
                continue;
            }
            auto y_rem = -ybeam;
            auto pT_i = std::sqrt(p_i[1]*p_i[1] + p_i[2]*p_i[2]);
            auto virtuality = p_i[0]*p_i[0] - pT_i*pT_i - p_i[3]*p_i[3];
            if (virtuality >= mass_min*mass_min) {
                y_rem = 0.5*log((p_i[0] + p_i[3])/(p_i[0] - p_i[3]));
            } else {
                //std::cout<< " Warning: The beam remnant of target is space-like vector." << std::endl;
                auto pT_temp = std::min(pT_i, std::sqrt(p_i[0]*p_i[0]-mass_min*mass_min));
                if (pT_temp < pT_i) {
                    //std::cout<< " Warning: E < pT of remnant of projectile, change pT from " << pT_i;
                    //std::cout<< " to " << pT_temp << std::endl;
                    p_i[1] = p_i[1]*pT_temp/pT_i;
                    p_i[2] = p_i[2]*pT_temp/pT_i;
                }
                //std::cout<< " Now change beam pz from " <<p_i[3];
                pT_i   = pT_temp;
                auto ymax = acosh(p_i[0] / std::sqrt(pT_i*pT_i + mass_min*mass_min));
                y_rem = -std::min(ymax, ybeam);
                auto mT_square = p_i[0]*p_i[0]/cosh(y_rem)/cosh(y_rem);
                p_i[3] = std::sqrt(mT_square) * sinh(y_rem);
                //std::cout<< " to " << p_i[3] << std::endl;
            }
            MomentumVec proj_p_vec = {p_i[0], p_i[1], p_i[2], -p_i[3]};

            get_tau_form_and_moversigma(string_evolution_mode, std::abs(y_rem),
                                        tau_form, m_over_sigma, y_loss);
            bool has_baryon_left = itarg->is_remnant_carry_baryon_number();
            bool has_baryon_right = false;
            QCDString qcd_string(x_i, tau_form, itarg, itarg,
                                 proj_p_vec, p_i, m_over_sigma,
                                 has_baryon_right, has_baryon_left);
            qcd_string.set_has_remnant_left(true);
            qcd_string.evolve_QCD_string();
            qcd_string.set_final_baryon_rapidities(y_rem + y_loss, 0.);
            remnant_string_list_.push_back(qcd_string);
        }
    }
}


void Glauber::update_collision_schedule(
                                shared_ptr<CollisionEvent> event_happened) {
    auto proj = event_happened->get_proj_nucleon_ptr().lock();
    proj->increment_collided_times();
    for (auto &it: (*proj->get_collide_nucleon_list()))
        create_a_collision_event(proj, it.lock());
    auto targ = event_happened->get_targ_nucleon_ptr().lock();
    targ->increment_collided_times();
    for (auto &it: (*targ->get_collide_nucleon_list()))
        create_a_collision_event(it.lock(), targ);
}


void Glauber::computeCenterOfMass(real &x_o, real &y_o) {
    // compute the center of mass
    x_o = 0.;
    y_o = 0.;
    int n_strings = QCD_string_list.size();
    for (auto &it: QCD_string_list) {
        auto x_prod = it.get_x_production();
        x_o += x_prod[1];
        y_o += x_prod[2];
    }
    if (sample_valence_quark) {
        for (auto &it: remnant_string_list_) {
            auto x_prod = it.get_x_production();
            x_o += x_prod[1];
            y_o += x_prod[2];
        }
        n_strings += remnant_string_list_.size();
    }
    x_o /= n_strings;
    y_o /= n_strings;
}


void Glauber::prepare_output_QCD_strings() {
    // compute the center of mass
    real x_o = 0.;
    real y_o = 0.;
    computeCenterOfMass(x_o, y_o);
    fill_QCD_string_output(QCD_string_list, remnant_string_list_, x_o, y_o,
                           QCD_string_output_arr_);
}


void Glauber::fill_QCD_string_output(std::vector<QCDString> &strings,
                                     std::vector<QCDString> &remnants,
                                     const real x_o, const real y_o,
                                     std::vector<std::vector<real>> &out) const {
    out.clear();
    // prepare output strings
    for (auto &it: strings) {
        auto x_prod = it.get_x_production();
        auto x_left = it.get_targ()->get_x();
        auto x_right = it.get_proj()->get_x();
        if (sample_valence_quark) {
            auto xq_left = it.get_targ_q()->get_x();
            x_left[1]  += xq_left[1];
            x_left[2]  += xq_left[2];

            auto xq_right = it.get_proj_q()->get_x();
            x_right[1] += xq_right[1];
            x_right[2] += xq_right[2];
        }
        auto tau_0  = sqrt(x_prod[0]*x_prod[0] - x_prod[3]*x_prod[3]);
        auto etas_0 = 0.5*log((x_prod[0] + x_prod[3])/(x_prod[0] - x_prod[3]));

        real remnant_left  = 0.;
        if (it.get_has_remnant_left()) {
            remnant_left = 1.0;
        }

        real remnant_right = 0.;
        if (it.get_has_remnant_right()) {
            remnant_right = 1.0;
        }

        real baryon_fraction_left  = 0.;
        if (it.get_has_baryon_left()) {
            baryon_fraction_left = 1.0;
        }

        real baryon_fraction_right = 0.;
        if (it.get_has_baryon_right()) {
            baryon_fraction_right = 1.0;
        }

        auto mass = it.get_string_mass();
        std::vector<real> output_array = {
            mass, it.get_m_over_sigma(), it.get_tau_form(),
            tau_0, etas_0, x_prod[1] - x_o, x_prod[2] - y_o,
            x_left[1] - x_o, x_left[2] - y_o,
            x_right[1] - x_o, x_right[2] - y_o,
            it.get_eta_s_left(), it.get_eta_s_right(),
            it.get_y_f_left(), it.get_y_f_right(),
            remnant_left, remnant_right,
            it.get_y_i_left(), it.get_y_i_right(),
            it.get_eta_s_baryon_left(), it.get_eta_s_baryon_right(),
            it.get_y_f_baryon_left(), it.get_y_f_baryon_right(),
            baryon_fraction_left, baryon_fraction_right,
            it.get_string_px(), it.get_string_py(),
        };
        out.push_back(output_array);
    }

    // output the beam remnant strings
    if (sample_valence_quark) {
        for (auto &it: remnants) {
            auto x_prod = it.get_x_production();
            auto x_left = it.get_targ()->get_x();
            auto x_right = it.get_proj()->get_x();
            auto tau_0  = sqrt(x_prod[0]*x_prod[0] - x_prod[3]*x_prod[3]);
            auto etas_0 = 0.5*log((x_prod[0] + x_prod[3])
                                  /(x_prod[0] - x_prod[3]));

            real remnant_left  = 0.;
            if (it.get_has_remnant_left()) {
                remnant_left = 1.0;
            }

            real remnant_right = 0.;
            if (it.get_has_remnant_right()) {
                remnant_right = 1.0;
            }

            real baryon_fraction_left  = 0.;
            if (it.get_has_baryon_left()) {
                baryon_fraction_left = 1.0;
            }

            real baryon_fraction_right = 0.;
            if (it.get_has_baryon_right()) {
                baryon_fraction_right = 1.0;
            }

            auto mass = it.get_string_mass();
            auto eta_s_center = (it.get_eta_s_left() + it.get_eta_s_right())/2.;
            auto eta_s_left = (  remnant_left*it.get_eta_s_left()
                               + (1. - remnant_left)
                                 *(std::max(
                                     (eta_s_center + it.get_eta_s_right())/2.,
                                     it.get_eta_s_right() - 1.0)));
            auto eta_s_right = (  remnant_right*it.get_eta_s_right()
                                + (1. - remnant_right)
                                  *(std::min(
                                     (eta_s_center + it.get_eta_s_left())/2.,
                                     it.get_eta_s_left() + 1.0)));
            std::vector<real> output_array = {
                mass, it.get_m_over_sigma(), it.get_tau_form(),
                tau_0, etas_0, x_prod[1] - x_o, x_prod[2] - y_o,
                x_left[1] - x_o, x_left[2] - y_o,
                x_right[1] - x_o, x_right[2] - y_o,
                eta_s_left, eta_s_right,
                remnant_left*it.get_y_f_left() + (1. - remnant_left)*eta_s_left,//13
                remnant_right*it.get_y_f_right() + (1. - remnant_right)*eta_s_right,//14
                remnant_left, remnant_right,//15 16 weight 0 or 1
                remnant_left*it.get_y_i_left() + (1. - remnant_left)*eta_s_left, //17 
                remnant_right*it.get_y_i_right() + (1. - remnant_right)*eta_s_right,//18
                remnant_left*it.get_eta_s_baryon_left(),
                remnant_right*it.get_eta_s_baryon_right(),
                remnant_left*it.get_y_f_baryon_left(),
                remnant_right*it.get_y_f_baryon_right(),
                baryon_fraction_left, baryon_fraction_right,
                it.get_string_px(), it.get_string_py(),
            };
            out.push_back(output_array);
        }
    }
}


namespace {

//! E and P_z of output rows (prepare_output_QCD_strings), as MUSIC computes
//! them (HydroSourceStrings: the strings' ends plus their remnants)
void output_energy(const std::vector<std::vector<real>> &out,
                   real &E, real &Pz) {
    E = 0.;
    Pz = 0.;
    for (const auto &r: out) {
        const real m = r[0];
        E  += (m*(cosh(r[17]) + cosh(r[18]) - cosh(r[13]) - cosh(r[14]))
               + r[15]*m*cosh(r[13]) + r[16]*m*cosh(r[14]));
        Pz += (m*(sinh(r[17]) + sinh(r[18]) - sinh(r[13]) - sinh(r[14]))
               + r[15]*m*sinh(r[13]) + r[16]*m*sinh(r[14]));
    }
}

//! the earliest start of the output rows' sources, by MUSIC's rule
//! (HydroSourceStrings: tau_0 + tau_form, or the end's tau if the string
//! lies on one side of its production point), which sets MUSIC's tau0
real tau_start(const real tau_0, const real tau_form, const real eta_0,
                const real eta_l, const real eta_r) {
    auto end_tau = [&](const real eta) {
        const real f1 = tau_0*tau_0 - tau_form*tau_form;
        const real f2 = tau_0*cosh(eta - eta_0);
        return(f2 + sqrt(f2*f2 - f1));
    };
    if (eta_l > eta_0) return(end_tau(eta_l));
    if (eta_r < eta_0) return(end_tau(eta_r));
    return(tau_0 + tau_form);
}

real min_tau_start(const std::vector<std::vector<real>> &out) {
    real tau_min = 1e30;
    for (const auto &r: out)
        tau_min = std::min(tau_min, tau_start(r[3], r[2], r[4], r[11], r[12]));
    return(tau_min);
}

//! A valence quark's strings on one side (the projectile's quarks are the
//! strings' right ends): its rapidity only drops (rises on the target side)
//! from string to string, so y0, its extreme y_i, gives its energy, and idx
//! runs from its first string to its last, which carries its remnant.
struct Chain {
    shared_ptr<Quark> q;
    shared_ptr<Nucleon> n;
    real y0;
    real mass;
    std::vector<std::size_t> idx;
    real energy() const {return(mass*cosh(y0));}
};

//! the chains of nucleon n's quarks (of every nucleon on that side if n is
//! null), in the order of their first string in strings
std::vector<Chain> quark_chains(std::vector<QCDString> &strings,
                                shared_ptr<Nucleon> n, const bool proj_side) {
    std::vector<Chain> chains;
    for (std::size_t i = 0; i < strings.size(); i++) {
        auto &s = strings[i];
        auto n_i = proj_side ? s.get_proj() : s.get_targ();
        if (n && n_i != n) continue;
        auto q = proj_side ? s.get_proj_q() : s.get_targ_q();
        const real y_i = proj_side ? s.get_y_i_right() : s.get_y_i_left();
        auto it = std::find_if(chains.begin(), chains.end(),
                               [&q](const Chain &c) {return c.q == q;});
        if (it == chains.end()) {
            chains.push_back({q, n_i, y_i, s.get_string_mass(), {i}});
        } else {
            it->idx.push_back(i);
            it->y0 = proj_side ? std::max(it->y0, y_i) : std::min(it->y0, y_i);
        }
    }
    for (auto &c: chains) {
        std::stable_sort(c.idx.begin(), c.idx.end(),
                         [&](std::size_t a, std::size_t b) {
                             return proj_side
                                 ? strings[a].get_y_i_right() > strings[b].get_y_i_right()
                                 : strings[a].get_y_i_left() < strings[b].get_y_i_left();});
    }
    return(chains);
}

//! one of the candidate chains, by a hash of their rapidities (no random
//! number: the event's stream is untouched)
std::size_t hash_pick(const std::vector<Chain> &chains,
                      const std::vector<std::size_t> &cand) {
    std::uint64_t h = 1469598103934665603ULL;
    for (auto k: cand) {
        std::uint64_t w;
        std::memcpy(&w, &chains[k].y0, sizeof(w));
        h = (h ^ w)*1099511628211ULL;
    }
    return(cand[h % cand.size()]);
}

}  // namespace


void Glauber::shift_string_end(QCDString &s, const bool proj_side,
                               const real dy, std::vector<real> &info) const {
    // evolve_QCD_string must give the string as produced
    QCDString check = s;
    check.evolve_QCD_string();
    if (std::abs(check.get_y_f_left()    - s.get_y_f_left())    > 1e-9
     || std::abs(check.get_y_f_right()   - s.get_y_f_right())   > 1e-9
     || std::abs(check.get_eta_s_left()  - s.get_eta_s_left())  > 1e-9
     || std::abs(check.get_eta_s_right() - s.get_eta_s_right()) > 1e-9)
        info[12] += 1;
    const real tau_form = s.get_tau_form();
    const real yb_l = s.get_y_f_baryon_left();
    const real yb_r = s.get_y_f_baryon_right();
    if (proj_side) {
        s.set_initial_rapidities(s.get_y_i_left(), s.get_y_i_right() + dy);
    } else {
        s.set_initial_rapidities(s.get_y_i_left() + dy, s.get_y_i_right());
    }
    s.evolve_QCD_string();
    if (s.get_tau_form() != tau_form) info[11] += 1;
    if (!parameter_list.get_baryon_junctions()) {
        s.set_final_baryon_rapidities(s.get_y_f_left(), s.get_y_f_right());
    } else {
        s.set_final_baryon_rapidities(yb_l + (proj_side ? 0. : dy),
                                      yb_r + (proj_side ? dy : 0.));
    }
    info[10] += 1;
}


real Glauber::max_quark_shift(std::vector<QCDString> &strings,
                              const std::vector<std::size_t> &idx,
                              const bool proj_side, const real dy_full,
                              const real tau_min) const {
    // the part of dy_full (from 0) by which a quark's strings can move without
    // a capped rapidity loss (which shortens tau_form), an end crossing the
    // production point, or a string starting before tau_min, the event's
    // earliest start: MUSIC's tau0 stays the event's
    auto ok = [&](const real dy) {
        for (auto i: idx) {
            QCDString s = strings[i];
            const real tau_form = s.get_tau_form();
            const auto x = s.get_x_production();
            const real eta_0 = 0.5*log((x[0] + x[3])/(x[0] - x[3]));
            const bool around = (s.get_eta_s_left() <= eta_0
                                 && eta_0 <= s.get_eta_s_right());
            if (proj_side) {
                s.set_initial_rapidities(s.get_y_i_left(),
                                         s.get_y_i_right() + dy);
            } else {
                s.set_initial_rapidities(s.get_y_i_left() + dy,
                                         s.get_y_i_right());
            }
            s.evolve_QCD_string();
            if (s.get_tau_form() != tau_form) return(false);
            if (around && !(s.get_eta_s_left() <= eta_0
                            && eta_0 <= s.get_eta_s_right()))
                return(false);
            const real tau_0 = sqrt(x[0]*x[0] - x[3]*x[3]);
            if (tau_start(tau_0, s.get_tau_form(), eta_0, s.get_eta_s_left(),
                          s.get_eta_s_right()) < tau_min)
                return(false);
        }
        return(true);
    };
    if (ok(dy_full)) return(dy_full);
    real lo = 0.;
    real hi = dy_full;
    for (int it = 0; it < 50; it++) {
        const real mid = 0.5*(lo + hi);
        if (ok(mid)) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return(ok(lo) ? lo : 0.);
}


void Glauber::take_energy_from_quarks(std::vector<QCDString> &strings,
                                      shared_ptr<Nucleon> n,
                                      const bool proj_side, const real E_take,
                                      const int i_side,
                                      std::vector<real> &info) const {
    const real tau_min = info[21];   // the event's earliest start
    auto chains = quark_chains(strings, n, proj_side);
    info[i_side] = E_take;
    if (chains.empty() || !(E_take > 0.)) return;

    // what each quark can give: at most 90 % of its energy above m_q, and no
    // more than keeps its strings' formation times (max_quark_shift)
    std::vector<real> E_max(chains.size());
    for (std::size_t k = 0; k < chains.size(); k++) {
        const auto &c = chains[k];
        const real sign = proj_side ? 1. : -1.;
        const real E_q = c.energy();
        const real dy_full = (sign*acosh((E_q - 0.9*(E_q - c.mass))/c.mass)
                              - c.y0);
        const real dy = max_quark_shift(strings, c.idx, proj_side, dy_full,
                                        tau_min);
        E_max[k] = E_q - c.mass*cosh(c.y0 + dy);
    }
    // one quark that can give E_take alone; else several, the most
    // energetic first
    std::vector<std::pair<std::size_t, real>> takes;
    std::vector<std::size_t> cand;
    for (std::size_t k = 0; k < chains.size(); k++)
        if (E_max[k] >= E_take) cand.push_back(k);
    if (!cand.empty()) {
        takes.push_back({hash_pick(chains, cand), E_take});
    } else {
        std::vector<std::size_t> order(chains.size());
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(),
                         [&](std::size_t a, std::size_t b) {
                             return chains[a].energy() > chains[b].energy();});
        real left = E_take;
        for (auto k: order) {
            if (!(left > 0.)) break;
            const real t = std::min(left, E_max[k]);
            if (!(t > 0.)) continue;
            takes.push_back({k, t});
            left -= t;
        }
    }

    real taken = 0.;
    for (std::size_t j = 0; j < takes.size(); j++) {
        auto &c = chains[takes[j].first];
        const real E_q = c.energy();
        const real y0_new = ((proj_side ? 1. : -1.)
                             *acosh((E_q - takes[j].second)/c.mass));
        const real dy = y0_new - c.y0;
        if (j == 0) {
            info[6 + i_side] = dy;
            info[8 + i_side] = E_q;
        }
        for (auto i: c.idx) shift_string_end(strings[i], proj_side, dy, info);
        taken += takes[j].second;
    }
    info[2 + i_side] = taken;
    info[4 + i_side] = takes.size();
}


real Glauber::remnant_rapidity(MomentumVec &p_i, const bool proj_side) const {
    // as produce_remnant_strings: the remnant's rapidity, a space-like
    // remnant put on the light cone (p_i changed accordingly)
    const real mass_min = 1.e-3;
    real y_rem = proj_side ? ybeam : -ybeam;
    auto pT_i = std::sqrt(p_i[1]*p_i[1] + p_i[2]*p_i[2]);
    auto virtuality = p_i[0]*p_i[0] - pT_i*pT_i - p_i[3]*p_i[3];
    if (virtuality >= mass_min*mass_min) {
        y_rem = 0.5*log((p_i[0] + p_i[3])/(p_i[0] - p_i[3]));
    } else {
        auto pT_temp = std::min(pT_i,
                                std::sqrt(p_i[0]*p_i[0] - mass_min*mass_min));
        if (pT_temp < pT_i) {
            p_i[1] = p_i[1]*pT_temp/pT_i;
            p_i[2] = p_i[2]*pT_temp/pT_i;
        }
        pT_i = pT_temp;
        auto ymax = acosh(p_i[0]/std::sqrt(pT_i*pT_i + mass_min*mass_min));
        y_rem = (proj_side ? 1. : -1.)*std::min(ymax, ybeam);
        auto mT_square = p_i[0]*p_i[0]/cosh(y_rem)/cosh(y_rem);
        p_i[3] = std::sqrt(mT_square)*sinh(y_rem);
    }
    return(y_rem);
}


int Glauber::set_remnant_string(std::vector<QCDString> &remnants,
                                shared_ptr<Nucleon> n, const bool proj_side,
                                const real E, std::vector<real> &info) const {
    // nucleon n's remnant string set to deposit E as MUSIC counts a remnant
    // (mass cosh y): an existing one keeps its rapidity and only its mass
    // scales, so its space-time evolution and start stay the event's (1; 3:
    // removed, E <= 1e-3 GeV); a nucleon without one gets one at its
    // remnant's freeze-out point with the settings of the first remnant
    // string of its side (2); 0: nothing to do
    const real mass_min = 1.e-3;
    auto it = std::find_if(remnants.begin(), remnants.end(),
                           [&](QCDString &r) {
                               return(r.get_proj() == n && r.get_targ() == n
                                      && (proj_side ? r.get_has_remnant_right()
                                                    : r.get_has_remnant_left()));});
    if (it != remnants.end()) {
        if (E <= mass_min) {
            remnants.erase(it);
            return(3);
        }
        const real y = proj_side ? it->get_y_i_right() : it->get_y_i_left();
        it->set_string_mass(E/cosh(y));
        return(1);
    }
    if (E <= mass_min) return(0);
    const QCDString *like = nullptr;
    for (auto &r: remnants)
        if (proj_side ? r.get_has_remnant_right() : r.get_has_remnant_left()) {
            like = &r;
            break;
        }
    if (!like) return(0);
    // beam-like, as the remnants are: on the light cone at the beam rapidity
    MomentumVec p_i = {E, 0., 0., 0.};
    const real y = proj_side ? ybeam : -ybeam;
    p_i[3] = E*tanh(y);
    MomentumVec mirror = {p_i[0], 0., 0., -p_i[3]};
    const bool baryon = n->is_remnant_carry_baryon_number();
    QCDString r = proj_side
        ? QCDString(n->get_remnant_x_frez(), like->get_tau_form(), n, n,
                    p_i, mirror, like->get_m_over_sigma(), baryon, false)
        : QCDString(n->get_remnant_x_frez(), like->get_tau_form(), n, n,
                    mirror, p_i, like->get_m_over_sigma(), false, baryon);
    if (proj_side) {
        r.set_has_remnant_right(true);
    } else {
        r.set_has_remnant_left(true);
    }
    r.evolve_QCD_string();
    if (r.get_tau_form() != like->get_tau_form()) info[11] += 1;
    if (proj_side) {
        r.set_final_baryon_rapidities(0., r.get_y_f_right());
    } else {
        r.set_final_baryon_rapidities(r.get_y_f_left(), 0.);
    }
    remnants.push_back(r);
    return(2);
}


void Glauber::remove_hard_quarks(std::vector<QCDString> &strings,
                                 std::vector<QCDString> &remnants,
                                 shared_ptr<Nucleon> proj_n,
                                 shared_ptr<Nucleon> targ_n,
                                 const MomentumVec &p_proj,
                                 const MomentumVec &p_targ,
                                 std::vector<real> &info) const {
    // 3dMCGlauber's ISR picture (Pick_and_subtract_hard_parton_momentum)
    // without iColored: the hard pair makes no string with each other, the
    // two hard quarks make no strings at all, and the energy they and the
    // pair's strings would have carried goes back to the nucleons' remnants,
    // from which the hard process's momentum is taken
    std::vector<char> rm(strings.size(), 0);
    shared_ptr<Quark> hard_q[2];
    for (int side = 0; side < 2; side++) {
        const bool proj_side = side == 0;
        auto chains = quark_chains(strings, proj_side ? proj_n : targ_n,
                                   proj_side);
        info[side] = (proj_side ? p_proj : p_targ)[0];
        if (chains.empty()) continue;
        std::vector<std::size_t> cand;
        for (std::size_t k = 0; k < chains.size(); k++)
            if (chains[k].energy() >= info[side]) cand.push_back(k);
        std::size_t k;
        if (!cand.empty()) {
            k = hash_pick(chains, cand);
        } else {
            k = 0;
            for (std::size_t j = 1; j < chains.size(); j++)
                if (chains[j].energy() > chains[k].energy()) k = j;
        }
        hard_q[side] = chains[k].q;
        info[8 + side] = chains[k].energy();
        for (auto i: chains[k].idx) rm[i] = 1;
        info[4 + side] = chains[k].idx.size();
    }
    for (std::size_t i = 0; i < strings.size(); i++)
        if (strings[i].get_proj() == proj_n && strings[i].get_targ() == targ_n)
            rm[i] = 1;

    // every other quark that loses a string keeps its energy: its later
    // strings start higher by the rapidity loss of the removed ones, and its
    // remnant moves to its last string kept; a quark left without strings
    // gives its energy back to its nucleon's remnant
    struct Back {
        shared_ptr<Nucleon> n;
        bool proj_side;
        MomentumVec p;
    };
    std::vector<Back> back;     // in the order met: the output's order fixed
    auto back_of = [&](shared_ptr<Nucleon> n, bool proj_side) -> Back & {
        auto it = std::find_if(back.begin(), back.end(),
                               [&n](const Back &b) {return b.n == n;});
        if (it != back.end()) return(*it);
        back.push_back({n, proj_side, {0., 0., 0., 0.}});
        return(back.back());
    };
    auto give_back = [&](const Chain &c, bool proj_side) {
        auto &b = back_of(c.n, proj_side);
        b.p[0] += c.energy();
        b.p[3] += c.mass*sinh(c.y0);
    };
    std::vector<std::tuple<std::size_t, bool, real>> shifts;
    std::vector<std::pair<std::size_t, bool>> new_remnant_end;
    for (int side = 0; side < 2; side++) {
        const bool proj_side = side == 0;
        for (auto &c: quark_chains(strings, nullptr, proj_side)) {
            bool touched = false;
            for (auto i: c.idx) touched = touched || rm[i];
            if (!touched) continue;
            if (c.q == hard_q[side]) {
                give_back(c, proj_side);
                continue;
            }
            real shift = 0.;
            std::size_t last_kept = strings.size();
            for (auto i: c.idx) {
                auto &s = strings[i];
                if (rm[i]) {
                    shift += (proj_side ? s.get_y_i_right() - s.get_y_f_right()
                                        : s.get_y_i_left() - s.get_y_f_left());
                } else {
                    if (shift != 0.) shifts.emplace_back(i, proj_side, shift);
                    last_kept = i;
                }
            }
            if (!rm[c.idx.back()]) continue;
            if (last_kept < strings.size()) {
                new_remnant_end.emplace_back(last_kept, proj_side);
            } else {
                give_back(c, proj_side);
            }
        }
    }
    for (auto &t: shifts)
        shift_string_end(strings[std::get<0>(t)], std::get<1>(t),
                         std::get<2>(t), info);
    for (auto &e: new_remnant_end) {
        if (e.second) {
            strings[e.first].set_has_remnant_right(true);
        } else {
            strings[e.first].set_has_remnant_left(true);
        }
    }

    // the remnants: what came back, the hard pair's minus the hard process
    // in MUSIC's measure of a remnant row (mass cosh y_i), so the strings
    // deposit exactly what is taken less (P_z follows the remnants'
    // rapidities)
    back_of(proj_n, true);
    back_of(targ_n, false);
    for (auto &b: back) {
        real E = 0.;
        auto it = std::find_if(remnants.begin(), remnants.end(),
                               [&](QCDString &r) {
                                   return(r.get_proj() == b.n && r.get_targ() == b.n
                                          && (b.proj_side ? r.get_has_remnant_right()
                                                          : r.get_has_remnant_left()));});
        if (it != remnants.end()) {
            const real y = b.proj_side ? it->get_y_i_right() : it->get_y_i_left();
            E = it->get_string_mass()*cosh(y);
        }
        const real E_before = E;
        E += b.p[0];
        int side = -1;
        if (b.n == proj_n) side = 0;
        if (b.n == targ_n) side = 1;
        if (side >= 0) {
            const MomentumVec &p_take = side == 0 ? p_proj : p_targ;
            const real available = std::max(0., E - 1.e-3);
            const real taken = std::min(p_take[0], available);
            const real f = p_take[0] > 0. ? taken/p_take[0] : 0.;
            E -= f*p_take[0];
            info[2 + side] = taken;
            info[6 + side] = E - E_before;
        }
        if (set_remnant_string(remnants, b.n, b.proj_side, E, info) != 0)
            info[19] += 1;
    }

    // a removed string stays as a placeholder without energy (no rapidity
    // loss, no remnant, no baryon number): it deposits nothing but keeps its
    // start, so the strings' earliest start, and with it MUSIC's tau0, stays
    // the event's
    for (std::size_t i = 0; i < strings.size(); i++) {
        if (!rm[i]) continue;
        auto &s = strings[i];
        s.set_final_rapidities(s.get_y_i_left(), s.get_y_i_right());
        s.set_has_remnant_left(false);
        s.set_has_remnant_right(false);
        if (s.get_has_baryon_left()) info[20] += 1;
        if (s.get_has_baryon_right()) info[20] += 1;
        s.set_has_baryon_left(false);
        s.set_has_baryon_right(false);
        info[18] += 1;
    }
}


std::vector<std::vector<real>>
Glauber::get_QCD_strings_output_list_hard_subtracted(
        const int collision_idx, const MomentumVec &p_proj,
        const MomentumVec &p_targ, const int mode, std::vector<real> &info) {
    info.assign(HARD_SUBTRACTION_INFO, 0.);
    if (!sample_valence_quark) {
        throw std::runtime_error(
            "Glauber: the hard subtraction needs valence quarks");
    }
    if (mode != 1 && mode != 2) {
        throw std::invalid_argument(
            "Glauber: unknown hard subtraction mode " + std::to_string(mode));
    }
    if (collision_idx < 0
        || collision_idx >= static_cast<int>(collision_schedule_list_.size())) {
        throw std::out_of_range(
            "Glauber: no binary collision " + std::to_string(collision_idx));
    }
    auto proj_n = collision_schedule_list_[collision_idx]
                                            .get_proj_nucleon_ptr().lock();
    auto targ_n = collision_schedule_list_[collision_idx]
                                            .get_targ_nucleon_ptr().lock();

    // the event's centre of mass: every other string stays where it is
    real x_o = 0.;
    real y_o = 0.;
    computeCenterOfMass(x_o, y_o);
    std::vector<QCDString> strings(QCD_string_list);
    std::vector<QCDString> remnants(remnant_string_list_);
    for (auto &s: strings)
        if (s.get_proj() == proj_n && s.get_targ() == targ_n) info[17] += 1;

    std::vector<std::vector<real>> out;
    fill_QCD_string_output(strings, remnants, x_o, y_o, out);
    output_energy(out, info[13], info[15]);
    info[21] = min_tau_start(out);
    if (mode == 1) {
        take_energy_from_quarks(strings, proj_n, true,  p_proj[0], 0, info);
        take_energy_from_quarks(strings, targ_n, false, p_targ[0], 1, info);
    } else {
        remove_hard_quarks(strings, remnants, proj_n, targ_n, p_proj, p_targ,
                           info);
    }
    fill_QCD_string_output(strings, remnants, x_o, y_o, out);
    output_energy(out, info[14], info[16]);
    info[22] = min_tau_start(out);
    return(out);
}


void Glauber::output_QCD_strings(std::string filename, const real Npart,
                                 const real Ncoll, const real Nstrings,
                                 const real b, const unsigned int seed) {
    if (QCD_string_output_arr_.size() == 0)
        prepare_output_QCD_strings();

    std::ofstream output(filename.c_str());
    real total_energy = Npart*parameter_list.get_roots()/2.;
    real net_Pz = ((projectile->get_number_of_wounded_nucleons()
                    - target->get_number_of_wounded_nucleons())
                   *parameter_list.get_roots()/2.);
    output << "# b = " << b << " fm " << "Npart = " << Npart
           << " Ncoll = " << Ncoll << " Nstrings = " << Nstrings
           << " total_energy = " << total_energy << " GeV, "
           << "net_Pz = " << net_Pz << " GeV, "
           << "seed = " << seed << endl;

    output << "# mass[GeV]  m_over_sigma[fm]  tau_form[fm]  tau_0[fm]  eta_s_0  "
           << "x_perp[fm]  y_perp[fm]  x_l[fm]  y_l[fm]  x_r[fm]  y_r[fm]  "
           << "eta_s_left  eta_s_right  y_l  y_r  remnant_l  remnant_r "
           << "y_l_i  y_r_i "
           << "eta_s_baryon_left  eta_s_baryon_right  y_l_baryon  y_r_baryon  "
           << "baryon_fraction_l  baryon_fraction_r"
           << endl;

    for (auto &string_i : QCD_string_output_arr_) {
        output << std::scientific << std::setprecision(8);
        for (auto &ival : string_i) {
            output << std::setw(15) << ival << "  ";
        }
        output << endl;
    }
    output.close();
}


void Glauber::prepareParticipantList() {
    participantList_.clear();
    // compute the center of mass
    real x_o = 0.;
    real y_o = 0.;
    computeCenterOfMass(x_o, y_o);
    auto proj_nucleon_list = projectile->get_nucleon_list();
    int dir = 1;
    for (auto &iproj: (*proj_nucleon_list)) {
        if (iproj->is_wounded()) {
            auto proj_x = iproj->get_x();
            proj_x[1] -= x_o;
            proj_x[2] -= y_o;
            std::vector<real> part_i;
            part_i.push_back(proj_x[0]);
            part_i.push_back(proj_x[1]);
            part_i.push_back(proj_x[2]);
            part_i.push_back(proj_x[3]);
            part_i.push_back(dir);
            part_i.push_back(iproj->get_electric_charge());
            participantList_.push_back(part_i);
        }
    }
    dir = -1;
    auto targ_nucleon_list = target->get_nucleon_list();
    for (auto &itarg: (*targ_nucleon_list)) {
        if (itarg->is_wounded()) {
            auto targ_x = itarg->get_x();
            targ_x[1] -= x_o;
            targ_x[2] -= y_o;
            std::vector<real> part_i;
            part_i.push_back(targ_x[0]);
            part_i.push_back(targ_x[1]);
            part_i.push_back(targ_x[2]);
            part_i.push_back(targ_x[3]);
            part_i.push_back(dir);
            part_i.push_back(itarg->get_electric_charge());
            participantList_.push_back(part_i);
        }
    }
}


void Glauber::outputParticipants(std::string filename) {
    prepareParticipantList();

    std::ofstream output(filename.c_str());
    output << "# t[fm]  x[fm]  y[fm]  z[fm]  dir  e"
           << endl;
    for (auto &ipart: participantList_) {
        output << std::scientific << std::setprecision(6);
        for (auto &x_i: ipart) {
            output << std::setw(10) << x_i << "  ";
        }
        output << endl;
    }
    output.close();
}


void Glauber::output_spectators(std::string filename) {
    // compute the center of mass
    real x_o = 0.;
    real y_o = 0.;
    computeCenterOfMass(x_o, y_o);

    std::ofstream output(filename.c_str());
    output << "# t[fm]  x[fm]  y[fm]  z[fm]  m[GeV]  px[GeV]  py[GeV]  y  e"
           << endl;
    auto proj_nucleon_list = projectile->get_nucleon_list();
    for (auto &iproj: (*proj_nucleon_list)) {
        if (!iproj->is_wounded()) {
            output << std::scientific << std::setprecision(6);
            auto proj_x = iproj->get_x();
            auto proj_p = iproj->get_p();
            auto mass = iproj->get_mass();

            // compute the when and where the spectator nucleon enters the
            // light cone
            real vz = proj_p[3]/proj_p[0];
            proj_x[3] = vz/(1. + vz)*(proj_x[3]/vz - proj_x[0]);
            proj_x[0] = -proj_x[3];
            proj_x[1] -= x_o;
            proj_x[2] -= y_o;

            // output spectator's position and momentum
            for (const auto &x_i : proj_x)
                output << std::setw(10) << x_i << "  ";

            auto fermiMomentum = iproj->get_fermi_momentum();
            for (int i = 1; i < 4; i++)
                proj_p[i] += fermiMomentum[i];
            auto rap = asinh(proj_p[3]/(sqrt(mass*mass + proj_p[1]*proj_p[1]
                                             + proj_p[2]*proj_p[2])));
            output << std::setw(10) << mass << "  "
                   << std::setw(10) << proj_p[1] << "  "
                   << std::setw(10) << proj_p[2] << "  "
                   << std::setw(10) << rap << "  "
                   << iproj->get_electric_charge() << endl;
        }
    }
    auto targ_nucleon_list = target->get_nucleon_list();
    for (auto &itarg: (*targ_nucleon_list)) {
        if (!itarg->is_wounded()) {
            output << std::scientific << std::setprecision(6);
            auto targ_x = itarg->get_x();
            auto targ_p = itarg->get_p();
            auto mass = itarg->get_mass();

            // compute the when and where the spectator nucleon enters the
            // light cone
            real vz = targ_p[3]/targ_p[0];
            targ_x[3] = vz/(1. - vz)*(targ_x[3]/vz - targ_x[0]);
            targ_x[0] = targ_x[3];
            targ_x[1] -= x_o;
            targ_x[2] -= y_o;

            // output spectator's position and momentum
            for (const auto &x_i : targ_x)
                output << std::setw(10) << x_i << "  ";

            auto fermiMomentum = itarg->get_fermi_momentum();
            for (int i = 1; i < 4; i++)
                targ_p[i] += fermiMomentum[i];
            auto rap = asinh(targ_p[3]/(sqrt(mass*mass + targ_p[1]*targ_p[1]
                                             + targ_p[2]*targ_p[2])));
            output << std::setw(10) << mass << "  "
                   << std::setw(10) << targ_p[1] << "  "
                   << std::setw(10) << targ_p[2] << "  "
                   << std::setw(10) << rap << "  "
                   << itarg->get_electric_charge() << endl;
        }
    }
    output.close();
}


real Glauber::sample_rapidity_loss_shell(real y_init) const {
    real y_loss = 0.0;
    if (parameter_list.get_rapidity_loss_method() == 1) {
        y_loss = sample_rapidity_loss_from_the_LEXUS_model(y_init);
    } else if (parameter_list.get_rapidity_loss_method() == 2) {
        y_loss = sample_rapidity_loss_from_parametrization(y_init);
    } else if (parameter_list.get_rapidity_loss_method() > 2) {
        y_loss = sample_rapidity_loss_from_parametrization_with_fluct(y_init);
    }
    return(y_loss);
}


real Glauber::sample_rapidity_loss_from_the_LEXUS_model(
                                            const real y_init) const {
    const real shape_coeff = 1.0;
    real sinh_y_lrf = sinh(shape_coeff*y_init);
    real arcsinh_factor = (ran_gen_ptr_->rand_uniform()
                           *(sinh(2.*shape_coeff*y_init) - sinh_y_lrf)
                           + sinh_y_lrf);
    real y_loss = 2.*y_init - asinh(arcsinh_factor)/shape_coeff;
    return(y_loss);
}


real Glauber::sample_rapidity_loss_from_parametrization(
                                                const real y_init) const {
    auto y_loss = (
        yloss_param_slope*pow(pow(y_init, yloss_param_a)*tanh(y_init),
                              yloss_param_b));
    return(y_loss);
}


real Glauber::sample_rapidity_loss_from_piecewise_parametrization(
                                                const real y_init) const {
    auto y_loss1 = parameter_list.getParam("ylossParam4At2", 1.60);
    auto y_loss2 = parameter_list.getParam("ylossParam4At4", 2.15);
    auto y_loss3 = parameter_list.getParam("ylossParam4At6", 2.45);
    auto y_loss4 = parameter_list.getParam("ylossParam4At10", 2.95);

    real y_loss = 0.;
    if (y_init < 2) {
        y_loss = y_loss1/2.*y_init;
    } else if (y_init < 4) {
        y_loss = (y_loss2 - y_loss1)/2.*y_init + (2.*y_loss1 - y_loss2);
    } else if (y_init < 6) {
        y_loss = (y_loss3 - y_loss2)/2.*y_init + (3.*y_loss2 - 2.*y_loss3);
    } else {
        y_loss = (y_loss4 - y_loss3)/4.*y_init + (2.5*y_loss3 - 1.5*y_loss4);
    }
    return(y_loss);
}


real Glauber::sample_rapidity_loss_from_parametrization_with_fluct(
                                                const real y_init) const {
    real y_mean = y_init/2.;
    real var = 0.5;
    if (parameter_list.get_rapidity_loss_method() == 3) {
        y_mean = sample_rapidity_loss_from_parametrization(y_init);
        auto y_rhic = 5.5;
        auto y_lhc = 9.0;
        auto var_rhic = parameter_list.get_yloss_param_fluct_var_RHIC();
        auto var_lhc = parameter_list.get_yloss_param_fluct_var_LHC();
        var = var_rhic;
        if (y_init > y_lhc) {
            var = var_lhc;
        } else if (y_init > y_rhic) {
            var = (var_rhic
                   + (y_init - y_rhic)/(y_lhc - y_rhic)*(var_lhc - var_rhic));
        }
    } else if (parameter_list.get_rapidity_loss_method() == 4) {
        y_mean = sample_rapidity_loss_from_piecewise_parametrization(y_init);
        var = parameter_list.getParam("ylossParam4var", 0.5);
    }

    // sample logit distribution for y_loss given the mean and variance
    auto random_x = ran_gen_ptr_->rand_normal(0., var);
    real logit_rand = 1./(1. + exp(-random_x));
    auto y_loss = y_mean;
    if (std::abs(2.*y_mean - y_init) > 1e-15) {
        real aa = (2.*y_mean - y_init)/(2.*y_mean*y_init*(y_init - y_mean));
        real bb = ((y_init*y_init - 2.*y_mean*y_mean)
                   /(2.*y_init*y_mean*(y_init - y_mean)));
        real cc = - logit_rand;
        y_loss = (-bb + sqrt(bb*bb - 4.*aa*cc))/(2.*aa);
    } else {
        y_loss = logit_rand*y_init;
    }
    return(y_loss);
}


// sample y from exp[(y - (0.5*(yt + yp)))/2]/(4.*Sinh[0.25*yp - 0.25*yt]),
// the new rapidity of the baryon number from the right moving particle
// after the collision in the lab frame
real Glauber::sample_junction_rapidity_right(real y_left, real y_right) const {
    real y = -2.*(-0.25*y_right - 0.25*y_left
                  - 1.*log(2.*(ran_gen_ptr_->rand_uniform()
                               + 0.5*exp(-0.25*y_right + 0.25*y_left)
                                 /sinh(0.25*y_right - 0.25*y_left))
                           *sinh(0.25*y_right - 0.25*y_left)));
    return(y);
}

// sample y from exp[-(y - (0.5*(yt + yp)))/2]/(4.*Sinh[0.25*yp - 0.25*yt]),
// the new rapidity of the baryon number from the left moving particle
// after the collision in the lab frame
real Glauber::sample_junction_rapidity_left(real y_left, real y_right) const {
    real y = 2.*(0.25*y_right + 0.25*y_left
                 - log(2.*(ran_gen_ptr_->rand_uniform()
                           + 0.5*exp(-0.25*y_right + 0.25*y_left)
                             /sinh(0.25*y_right - 0.25*y_left))
                       *sinh( 0.25 * y_right - 0.25 * y_left)));
    return(y);
}


// sample y from a uniformed 1/(yp - yt) distribution
real Glauber::sample_junction_rapidity_uniformed(real y_left, real y_right) const {
    real y = y_left + ran_gen_ptr_->rand_uniform()*(y_right - y_left);
    return(y);
}


// This function computes the sigeff(s) from the formula
// sigmaNN_in(s) = int d^2b [1 - exp(-sigeff(s)*Tpp(b))]
// Reads sigmaNN, returns guassian width
real Glauber::get_sig_eff(const real siginNN) {
    // rms-radius of a gaussian = rms-radius of a disc with radius R,
    // where 2*PI*(2R)^2=sigmaNN
    const real width = sqrt(0.1*siginNN/M_PI)/sqrt(8.);
    nucleon_width_ = width;
    const real sigin = siginNN*0.1;   // sigma_in(s) [mb --> fm^2]

    const int Nint = 100;    // # of integration points
    std::vector<real> b(Nint, 0.);
    std::vector<real> Tnn(Nint, 0.);
    const real Bmax = 10.0*width;
    const real db = Bmax/Nint;
    for (int i = 0; i < Nint; i++) {
        b[i] = (i + 0.5)*db;
        Tnn[i] = exp(-b[i]*b[i]/(4.*width*width))/(M_PI*(4.*width*width));
    }
    const real prefactor = 2.*M_PI*db;

    real sum, dN;

    real sigeff = 10.0;        // starting point of iteration [fm^2]
    real sigeff0;              // holds value from previous iteration step
    int tol = 0;
    do {                       // iterate ...
        sigeff0 = sigeff;
        sum = 0.0;
        dN  = 0.0;
        for (int ib = 0; ib < Nint; ib++) {  // integral d^2b from 0 to Bmax
            sum += prefactor*b[ib]*(1.0 - exp(-sigeff*Tnn[ib]));
            dN  += prefactor*b[ib]*Tnn[ib]*exp(-sigeff*Tnn[ib]);
        }
        sigeff -= (sum - sigin)/dN;
        tol++;
    } while(std::abs(sigeff - sigeff0) > 1e-4 && tol < 1000);
    // until sigeff has converged

    return(sigeff);
}

}

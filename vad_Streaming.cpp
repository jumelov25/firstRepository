/**
 * VOCALIZE - Speech and Language Technology Solutions
  *
 *      "VAD for Streaming operation" 
 *
 * @author VOCALIZE Team
 * @version 9.7
 *
 * (All Rights Reserved)
 */
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "voc_vad.h"
#include "sub_band_energy_level.h"

// ==============================
//    VOICE ACTIVITY DETECTOR   
// ==============================

static void voc_activity_detector_state_change(voc_activity_detector_t *detector, voc_detector_state_e state);
static double voc_activity_detector_level_calculate(voc_activity_detector_t *detector, short *buffer, int len);
static double voc_level_calculate(voc_activity_detector_t *detector, short *x, int lsize, int fs,
                                  int energy, double alpha, double beta);

/** Create activity detector */
voc_activity_detector_t *voc_activity_detector_create(int audio_rate)
{
    if (audio_rate != 22050 && audio_rate != 16000 && audio_rate != 8000) {
        return NULL;
    }

    voc_activity_detector_t *detector = new voc_activity_detector_t();
    detector->level_threshold = 0.1;
    detector->speech_timeout = 80;
    detector->silence_timeout = 300;
    detector->recognition_timeout = 9999999;
    detector->no_input_timeout = 15000;
    detector->noinput_duration = 0;
    detector->speech_duration = 0;
    detector->recognition_duration = 0;
    detector->silence_duration = 0;
    detector->total_duration = 0;
    detector->state = DETECTOR_STATE_INACTIVITY;

    // Dynamic threshold initializations
    detector->dynamic_threshold = false;
    detector->n_points = 0;
    detector->audio_rate = audio_rate;
    detector->level_floor = 0.02;
    detector->level_ceil = 0.98;
    detector->histogram_memory_length = 10000;
    detector->acc_time = 0.0;
    detector->activation_timer = 0.0;
    detector->step_factor = 1e-2;
    detector->mean_factor = 1.0;
    detector->alpha = 1.0;
    detector->beta = 0.0;

    // Memory
    detector->memory_len_ec_sub_band = 0;

    // 10 é o tamanho da média móvel,
    // inicializa a primeira métade com zeros
    for (int i = 0; i < 5; i++) {
        detector->memory_nvad[i] = 0.0;
    }
    detector->memory_len_nvad = 5;

    // 101 é o número de bins do histograma
    for (int i = 0; i < 101; i++) {
        detector->level_histogram[i] = 0;
    }

    return detector;
}


void voc_activity_detector_destroy(voc_activity_detector_t *detector)
{
    delete detector;
}


/** Reset activity detector */
void voc_activity_detector_reset(voc_activity_detector_t *detector)
{
    detector->noinput_duration = 0;
    detector->speech_duration = 0;
    detector->recognition_duration = 0;
    detector->silence_duration = 0;
    detector->total_duration = 0;
    detector->state = DETECTOR_STATE_INACTIVITY;
}


void voc_activity_detector_level_floor_set(voc_activity_detector_t *detector, double level)
{
    if (level >= 0.0 && level <= 1.0) {
        detector->level_floor = level;
    }
}


void voc_activity_detector_level_ceil_set(voc_activity_detector_t *detector, double level)
{
    if (level >= 0.0 && level <= 1.0) {
        detector->level_ceil = level;
    }
}


void voc_activity_detector_histogram_memory_length_set(voc_activity_detector_t *detector, int points)
{
    if (points > 0) {
        detector->histogram_memory_length = points;
    }
}

void voc_activity_detector_dynamic_threshold_set(voc_activity_detector_t *detector, bool use_dynamic_threshold)
{
    detector->dynamic_threshold = use_dynamic_threshold;
}


/** Set threshold of voice activity (silence) level */
void voc_activity_detector_level_set(voc_activity_detector_t *detector, double level_threshold)
{
    //printf("level %lu\n", (unsigned long) level_threshold);
    if (level_threshold >= 0.0 && level_threshold <= 1.0) {
        detector->level_threshold = level_threshold;
    }
}


/** Set the step factor for the dynamic threshold method */
void voc_activity_detector_step_factor_set(voc_activity_detector_t *detector, double factor) {
    if (factor >= 1e-4 && factor <= 1) {
        detector->step_factor = factor;
    }
}


/** */
void voc_activity_detector_target_factor_set(voc_activity_detector_t *detector, double factor) {
    if (factor >= 0.0 && factor <= 1.0) {
        detector->mean_factor = factor;
    }
}


/** */
void voc_activity_detector_alpha_set(voc_activity_detector_t *detector, double alpha) {
    detector->alpha = alpha;
}


/** */
void voc_activity_detector_beta_set(voc_activity_detector_t *detector, double beta) {
    detector->beta = beta;
}

/** */
void voc_activity_detector_activation_timer_set(voc_activity_detector_t *detector, double timer) {
    if (timer > 0) {
        detector->activation_timer = timer;
    }
}

/** Set noinput timeout */
void voc_activity_detector_no_input_timeout_set(voc_activity_detector_t *detector, size_t no_input_timeout)
{
    //printf("no_input_timeout %lu\n", (unsigned long) no_input_timeout);
    detector->no_input_timeout = no_input_timeout;
}


/** Set timeout required to trigger speech (transition from inactive to active state) */
void voc_activity_detector_speech_timeout_set(voc_activity_detector_t *detector, size_t speech_timeout)
{
    //printf("speech_timeout %lu\n", (unsigned long) speech_timeout);
    detector->speech_timeout = speech_timeout;
}


/** Set timeout required to trigger silence (transition from active to inactive state) */
void voc_activity_detector_silence_timeout_set(voc_activity_detector_t *detector, size_t silence_timeout)
{
    //printf("silence_timeout %lu\n", (unsigned long) silence_timeout);
    detector->silence_timeout = silence_timeout;
}


/** Set recognition timeout */
void voc_activity_detector_recognition_timeout_set(voc_activity_detector_t *detector, size_t recognition_timeout)
{
    //printf("recognition_timeout %lu\n", (unsigned long) recognition_timeout);
    detector->recognition_timeout = recognition_timeout;
}


static void voc_activity_detector_state_change(voc_activity_detector_t *detector, voc_detector_state_e state)
{
    detector->state = state;
}


void voc_activity_detector_reset_durations(voc_activity_detector_t *detector, int dur_number)
{
    if (dur_number == 0) { // noinput_duration
        detector->noinput_duration = 0;
    }
    if (dur_number == 1) { // speech_duration
        detector->speech_duration = 0;
    }
    if (dur_number == 2) { // recognition_duration
        detector->recognition_duration = 0;
    }
    if (dur_number == 3) { // silence_duration
        detector->silence_duration = 0;
    }
}


static double voc_activity_detector_level_calculate(voc_activity_detector_t *detector, short *buffer, int len)
{
    int sum = 0;
    int mean = 0;
    int count = len;
    short *cur = buffer;
    short *sample_idx = cur;
    short *end = cur + count;

    // Cálculo da energia
    for (; sample_idx < end; sample_idx++){
        mean += *sample_idx;
    }
    mean = mean / count;

    for(; cur < end; cur++){
        if(*cur < 0) {
            sum -= (*cur - mean);
        }
        else {
            sum += (*cur - mean);
        }
    }

    // Cálculo da torção
    double level = voc_level_calculate(detector, buffer, len,
        detector->audio_rate, sum / count, detector->alpha, detector->beta);

    return level;
}


/** Process current frame */
voc_detector_event_e voc_activity_detector_process(voc_activity_detector_t *detector, short *buffer, int len)
{
    voc_detector_event_e det_event = VOC_DETECTOR_EVENT_NONE;
    double level = 0;
    detector->total_duration += VOC_CODEC_FRAME_TIME_BASE;

    // check frame length
    if (len != VOC_CODEC_FRAME_TIME_BASE * detector->audio_rate / 1000) {
        det_event = VOC_DETECTOR_EVENT_ERROR;
        return det_event;
    }

    // calculate current activity level of processed frame
    level = voc_activity_detector_level_calculate(detector, buffer, len);

    if (detector->dynamic_threshold) {
        double step;
        double mean = 0.0;
        //double variance = 0.0;
        int ilevel = (int) (100 * level);

        detector->acc_time += VOC_CODEC_FRAME_TIME_BASE;

        if (detector->acc_time < detector->activation_timer) {
            level = 0.0;
        }

        if (ilevel >= 0 && ilevel <= 100) {
            detector->hist_track.push(histogram_point_t(detector->acc_time, ilevel));
            detector->level_histogram[ilevel]++;
            detector->n_points++;

            if (detector->n_points > detector->histogram_memory_length) {
 
                detector->level_histogram[detector->hist_track.top()._level]--;
                detector->hist_track.pop();
                detector->n_points--;
            }

            for (int i = 0; i < 101; i++) {
                mean += (double) (i * detector->level_histogram[i]);
              }
            mean /= detector->n_points;
 
            step = (detector->mean_factor * mean / 100.0 - detector->level_threshold) * detector->step_factor;

            if (detector->level_threshold + step > detector->level_floor &&
                detector->level_threshold + step < detector->level_ceil) {
                voc_activity_detector_level_set(detector, detector->level_threshold + step);
            }
        }

  
    }

    // Esta seção governa as transições entre estados da máquina de estados
    if (detector->state == DETECTOR_STATE_INACTIVITY) {
        if (level >= detector->level_threshold) {

            /* start to detect activity */
            detector->noinput_duration += VOC_CODEC_FRAME_TIME_BASE;
            detector->speech_duration += VOC_CODEC_FRAME_TIME_BASE;
            voc_activity_detector_state_change(detector,DETECTOR_STATE_ACTIVITY_TRANSITION);
        } else {
            detector->noinput_duration += VOC_CODEC_FRAME_TIME_BASE;
            if(detector->noinput_duration >= detector->no_input_timeout) {
                
                /* detected noinput */
                det_event = VOC_DETECTOR_EVENT_NOINPUT;
            }
        }
    } else if (detector->state == DETECTOR_STATE_ACTIVITY_TRANSITION) {
        if (level >= detector->level_threshold) {
            detector->noinput_duration += VOC_CODEC_FRAME_TIME_BASE;
            detector->speech_duration += VOC_CODEC_FRAME_TIME_BASE;
            detector->recognition_duration += VOC_CODEC_FRAME_TIME_BASE;
            if (detector->speech_duration >= detector->speech_timeout) {

                /* finally detected activity */
                detector->recognition_duration += VOC_CODEC_FRAME_TIME_BASE;
                det_event = VOC_DETECTOR_EVENT_ACTIVITY;
                voc_activity_detector_state_change(detector,DETECTOR_STATE_ACTIVITY);
                voc_activity_detector_reset_durations(detector, 0);
            }
        } else {

            /* fallback to inactivity */
            detector->noinput_duration += VOC_CODEC_FRAME_TIME_BASE;
            voc_activity_detector_state_change(detector,DETECTOR_STATE_INACTIVITY);
            voc_activity_detector_reset_durations(detector, 1); // zerando speech_duration
            voc_activity_detector_reset_durations(detector, 2); // zerando recognition_duration
        }
    } else if (detector->state == DETECTOR_STATE_ACTIVITY) {
        if (level >= detector->level_threshold) {
            detector->recognition_duration += VOC_CODEC_FRAME_TIME_BASE;
            voc_activity_detector_reset_durations(detector, 1); // zerando speech_duration
            if(detector->recognition_duration >= detector->recognition_timeout) {

                /* detected inactivity */
                det_event = VOC_DETECTOR_EVENT_RECOGNITION_TIMEOUT;
                voc_activity_detector_state_change(detector,DETECTOR_STATE_INACTIVITY);
            }
        } else {

            /* start to detect inactivity */
            detector->recognition_duration += VOC_CODEC_FRAME_TIME_BASE;
            voc_activity_detector_state_change(detector,DETECTOR_STATE_INACTIVITY_TRANSITION);
        }
    } else if (detector->state == DETECTOR_STATE_INACTIVITY_TRANSITION) {
        if (level >= detector->level_threshold) { // 4c

            /* fallback to activity */
            detector->recognition_duration += VOC_CODEC_FRAME_TIME_BASE;
            voc_activity_detector_state_change(detector,DETECTOR_STATE_ACTIVITY);
            voc_activity_detector_reset_durations(detector, 3); // zerando silence_duration
        } else {
            detector->recognition_duration += VOC_CODEC_FRAME_TIME_BASE;
            detector->silence_duration += VOC_CODEC_FRAME_TIME_BASE;
            if (detector->silence_duration >= detector->silence_timeout) {

                /* detected inactivity */
                det_event = VOC_DETECTOR_EVENT_INACTIVITY;
                voc_activity_detector_state_change(detector,DETECTOR_STATE_INACTIVITY);
            }
            if (detector->recognition_duration >= detector->recognition_timeout) {

                /* detected inactivity */
                det_event = VOC_DETECTOR_EVENT_RECOGNITION_TIMEOUT;
                voc_activity_detector_state_change(detector,DETECTOR_STATE_INACTIVITY);
            }
        }
    }    

    return det_event;
}

static double voc_level_calculate(voc_activity_detector_t *detector, short *x, int lsize,
                                  int fs, int energy, double alpha, double beta) {
    double mn2 = 0.0;
    double level = 0.0;
    double level_1 = 0.0;
    double level_2 = 0.0;
    double *ec_sub_band = (double*) malloc(sizeof(double) * 5);
    double *g2 = (double*) malloc(sizeof(double) * 5);
    double *frame_copy = (double *) malloc(sizeof(double) * lsize);

    for (int i = 0; i < lsize; i++) {
        frame_copy[i] = 1.0 * x[i];
    }

    // Cálculo da energia por sub-banda (tamanho da FFT: 256)
    sub_band_energy_level_vad(frame_copy, lsize, ec_sub_band, fs, 256);

    // Atualização da memória
    if (detector->memory_len_ec_sub_band < 4) {
        detector->memory_ec_sub_band[detector->memory_len_ec_sub_band][0] = ec_sub_band[0];
        detector->memory_ec_sub_band[detector->memory_len_ec_sub_band][1] = ec_sub_band[1];
        detector->memory_ec_sub_band[detector->memory_len_ec_sub_band][2] = ec_sub_band[2];
        detector->memory_ec_sub_band[detector->memory_len_ec_sub_band][3] = ec_sub_band[3];
        detector->memory_ec_sub_band[detector->memory_len_ec_sub_band][4] = ec_sub_band[4];

        detector->memory_len_ec_sub_band++;

        level = 0.0;
    } else {
        if (detector->memory_len_ec_sub_band == 4) {
            detector->memory_ec_sub_band[4][0] = ec_sub_band[0];
            detector->memory_ec_sub_band[4][1] = ec_sub_band[1];
            detector->memory_ec_sub_band[4][2] = ec_sub_band[2];
            detector->memory_ec_sub_band[4][3] = ec_sub_band[3];
            detector->memory_ec_sub_band[4][4] = ec_sub_band[4];

            detector->memory_len_ec_sub_band++;
        } else {
            for (int i = 1; i < 5; i++) {
                detector->memory_ec_sub_band[i - 1][0] = detector->memory_ec_sub_band[i][0];
                detector->memory_ec_sub_band[i - 1][1] = detector->memory_ec_sub_band[i][1];
                detector->memory_ec_sub_band[i - 1][2] = detector->memory_ec_sub_band[i][2];
                detector->memory_ec_sub_band[i - 1][3] = detector->memory_ec_sub_band[i][3];
                detector->memory_ec_sub_band[i - 1][4] = detector->memory_ec_sub_band[i][4];
            }

            detector->memory_ec_sub_band[4][0] = ec_sub_band[0];
            detector->memory_ec_sub_band[4][1] = ec_sub_band[1];
            detector->memory_ec_sub_band[4][2] = ec_sub_band[2];
            detector->memory_ec_sub_band[4][3] = ec_sub_band[3];
            detector->memory_ec_sub_band[4][4] = ec_sub_band[4];
        }

        // Cálculo do coef. de torção
        g2[0]  = (-2*detector->memory_ec_sub_band[0][0] -1*detector->memory_ec_sub_band[1][0] +0*detector->memory_ec_sub_band[2][0] +1*detector->memory_ec_sub_band[3][0] + 2*detector->memory_ec_sub_band[4][0]) / 10.0;
        g2[1]  = (-2*detector->memory_ec_sub_band[0][1] -1*detector->memory_ec_sub_band[1][1] +0*detector->memory_ec_sub_band[2][1] +1*detector->memory_ec_sub_band[3][1] + 2*detector->memory_ec_sub_band[4][1]) / 10.0;
        g2[2]  = (-2*detector->memory_ec_sub_band[0][2] -1*detector->memory_ec_sub_band[1][2] +0*detector->memory_ec_sub_band[2][2] +1*detector->memory_ec_sub_band[3][2] + 2*detector->memory_ec_sub_band[4][2]) / 10.0;
        g2[3]  = (-2*detector->memory_ec_sub_band[0][3] -1*detector->memory_ec_sub_band[1][3] +0*detector->memory_ec_sub_band[2][3] +1*detector->memory_ec_sub_band[3][3] + 2*detector->memory_ec_sub_band[4][3]) / 10.0;
        g2[4]  = (-2*detector->memory_ec_sub_band[0][4] -1*detector->memory_ec_sub_band[1][4] +0*detector->memory_ec_sub_band[2][4] +1*detector->memory_ec_sub_band[3][4] + 2*detector->memory_ec_sub_band[4][4]) / 10.0;

        for (int j = 0; j < 5; ++j){
            mn2 += g2[j]*g2[j];
        }

        // Normalização da torção
        level_1 = -mn2 / (6.0 * 256.0 * 16000.0 * 16000.0);
        level_1 = exp(level_1);
        level_1 = 2.0 / (1.0 + level_1) - 1.0;

        // Normalização da energia
        level_2 = -10 * energy / (VOC_CODEC_FRAME_TIME_BASE * 16000.0);
        level_2 = exp(level_2);
        level_2 = 2.0 / (1.0 + level_2) - 1.0;

        // Combinação: torção e energia
        level = -alpha * level_1 - beta * level_2;
        level = exp(level);
        level = 2.0 / (1.0 + level) - 1.0;
    }

    // Moving average
    double level_sum = 0.0;
    // 10 é o tamanho da média móvel
    if (detector->memory_len_nvad < 10) {
        detector->memory_nvad[detector->memory_len_nvad] = level;
        detector->memory_len_nvad++;
    } else {
        for (int i = 1; i < 10; i++) {
            detector->memory_nvad[i - 1] = detector->memory_nvad[i];
        }
        detector->memory_nvad[9] = level;
    }
    for (int i = 0; i < detector->memory_len_nvad; i++) {
        level_sum += detector->memory_nvad[i];
    }
    level = level_sum / detector->memory_len_nvad;

    free(frame_copy);
    free(ec_sub_band);
    free(g2);

    return level;
}

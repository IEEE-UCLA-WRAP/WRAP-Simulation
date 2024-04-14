/*
@param rrc_samps_in Pointer to buffer of RRC-filtered samples (after carrier recovery)
@param symbs_out Pointer to output buffer to place symbols from timing recovery
@param params Pointer to parameter struct which must have **TR_phase**, **TR_integrator**, and **sps** as properties
*/
int timing_pll(float* rrc_samps_in, float* symbs_out, params_r* params) {
    
    const float Kp = 2.5;
    const float Ki = 0.5;

    float error = 0;

    int i = 1;
    symbs_out[0] = rrc_samps_in[1*params->sps];

    while ((i < SYMBOL_BUFF) && (i*params->sps + params->TR_phase < ADC_BUF_LEN)) {
        symbs_out[i] = rrc_samps_in[i*params->sps + params->TR_phase];

        error = sign(symbs_out[i-1]) * symbs_out[i] - sign(symbs_out[i]) * symbs_out[i - 1];
        
        params->TR_integrator += error;
        params->TR_phase += Kp * error + Ki * params->TR_integrator;
        i++;
    }
    return i;
}

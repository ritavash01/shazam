#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/shm.h>
#include <sys/time.h>
#include "gmrt_newcorr.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <sys/time.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl.h> 
#define DasHeaderKey 2031

BeamHeaderType *dataHdr_FRB;

char *gtac_code; // we need to provide an argument for gtac code, but I am leaving it for now




namespace nb = nanobind;

nb::dict initialize_HDR_SHM_py() {
    nb::dict header_data; // Create a Python dictionary to store the data

    int id_FRB_hdr = shmget(DasHeaderKey, sizeof(BeamHeaderType), IPC_CREAT | 0666);
    if (id_FRB_hdr < 0) {
        perror("FRB shmget HEADER ID");
        return nb::dict(); // Return an empty dictionary on error
    }

    dataHdr_FRB = (BeamHeaderType *)shmat(id_FRB_hdr, 0, 0);
    if ((void *)dataHdr_FRB == (void *)-1) {
        perror("FRB shmat HEADER ID");
        return nb::dict();
    }

    ScanInfoType *Scan = nullptr;
    for (int i = 0; i < MAX_SCANS; i++) {
        Scan = &(dataHdr_FRB->ScanTab[i]);
        if (strcmp(Scan->proj.code, gtac_code) == 0) {
            break;
        }
    }

    if (Scan == nullptr) {
        fprintf(stderr, "GTAC code not found in ScanTab.\n");
        return nb::dict();
    }

    // Populate dictionary with header data
    header_data["GTAC_code"] = std::string(Scan->proj.code);
    header_data["Observer"] = std::string(Scan->proj.observer);
    header_data["Title"] = std::string(Scan->proj.title);
    header_data["Source"] = std::string(Scan->source.object);
    header_data["RA"] = Scan->source.ra_app;
    header_data["DEC"] = Scan->source.dec_app;
    header_data["Channels"] = dataHdr_FRB->corr.corrpar.channels;
    header_data["Bandwidth_MHz"] = dataHdr_FRB->corr.daspar.gsb_acq_bw / dataHdr_FRB->corr.daspar.gsb_final_bw;
    header_data["Frequency_Ch_0_Hz"] = Scan->source.freq[0];
    header_data["Channel_width_Hz"] = Scan->source.ch_width * -1 * Scan->source.net_sign[0];
    header_data["Sampling_time_uSec"] = dataHdr_FRB->corr.daspar.gsb_final_bw * dataHdr_FRB->BeamGenHdr.SampInterval /
                                        (dataHdr_FRB->corr.corrpar.clock / 1e6);

    int BeamID = dataHdr_FRB->BeamGenHdr.BeamHostID;
    header_data["Beam_ID"] = BeamID;
    header_data["Beam_mode"] = std::string(beam_type[dataHdr_FRB->BeamGenHdr.BeamType[BeamID]]);
    header_data["No_of_stokes"] = dataHdr_FRB->BeamGenHdr.NStokes[BeamID];
    header_data["Num_bits_per_sample"] = dataHdr_FRB->BeamGenHdr.OutputDataFormat;

    // Antenna data processing
    char ant_list[30][4] = {
        "C00", "C01", "C02", "C03", "C04", "C05", "C06", "C08", "C09", "C10",
        "C11", "C12", "C13", "C14", "E02", "E03", "E04", "E05", "E06", "S01",
        "S02", "S03", "S04", "S06", "W01", "W02", "W03", "W04", "W05", "W06"
    };

    unsigned int ref_ant_mask = 1;
    nb::list pol1_antennas;
    nb::list pol2_antennas;

    if (dataHdr_FRB->BeamGenHdr.BeamType[BeamID] == 1 || 
        dataHdr_FRB->BeamGenHdr.BeamType[BeamID] == 2 ||
        dataHdr_FRB->BeamGenHdr.BeamType[BeamID] == 4) {

        for (int ant_num = 0; ant_num < 30; ant_num++) {
            if ((ref_ant_mask << ant_num) & dataHdr_FRB->BeamGenHdr.GAC_maskP1[BeamID]) {
                pol1_antennas.append(std::string(ant_list[ant_num]));
            }
            if ((ref_ant_mask << ant_num) & dataHdr_FRB->BeamGenHdr.GAC_maskP2[BeamID]) {
                pol2_antennas.append(std::string(ant_list[ant_num]));
            }
        }

        header_data["Pol1_mask"] = std::to_string(dataHdr_FRB->BeamGenHdr.GAC_maskP1[BeamID]);
        header_data["Pol1_antennas"] = pol1_antennas;
        header_data["Pol2_mask"] = std::to_string(dataHdr_FRB->BeamGenHdr.GAC_maskP2[BeamID]);
        header_data["Pol2_antennas"] = pol2_antennas;

    } else if (dataHdr_FRB->BeamGenHdr.BeamType[BeamID] == 3) {
        for (int ant_num = 0; ant_num < 30; ant_num++) {
            if ((ref_ant_mask << ant_num) & dataHdr_FRB->BeamGenHdr.GAC_vlt_maskP1[BeamID]) {
                pol1_antennas.append(std::string(ant_list[ant_num]));
            }
            if ((ref_ant_mask << ant_num) & dataHdr_FRB->BeamGenHdr.GAC_vlt_maskP2[BeamID]) {
                pol2_antennas.append(std::string(ant_list[ant_num]));
            }
        }

        header_data["Pol1_mask"] = std::to_string(dataHdr_FRB->BeamGenHdr.GAC_vlt_maskP1[BeamID]);
        header_data["Pol1_antennas"] = pol1_antennas;
        header_data["Pol2_mask"] = std::to_string(dataHdr_FRB->BeamGenHdr.GAC_vlt_maskP2[BeamID]);
        header_data["Pol2_antennas"] = pol2_antennas;
    }

    return header_data;
}

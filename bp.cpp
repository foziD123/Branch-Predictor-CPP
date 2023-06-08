/* 046267 Computer Architecture - Winter 20/21 - HW #1                  */
/* This file should hold your implementation of the predictor simulator */

#include "bp_api.h"
#include <iostream>
#include <vector>
#include <math.h>


#define VALID_BIT_S 1
#define TARGET_PC 30
#define MAIN_TABLE_COL 4


using namespace std;
enum state {
    SNT = 0, WNT = 1, WT = 2, ST = 3
};
enum share {
    NOT_USING_SHARE, USING_SHARE_LSB, USING_SHARE_MID
};

typedef class BTB_t {
public:
    unsigned entries;           //the size of the BTB = the number of entries in the table
    unsigned BHR_size;          //size of the history register.
    unsigned tag_size;          //size of the tag
    unsigned start_state;       //the initial value of the state in the fsm table.
    bool isGlobalHist;          //global or local history register.
    bool isGlobalTable;         //global or local fsm
    uint32_t Shared;            //Lshare or Gshare (relevant for global table)


    vector<vector<uint32_t>> HistoryTable;       // matrix of unsigned int --> first col valid bit, sec col tag, 3rd col history, 4th col target PC
    vector<vector<uint32_t>> FSM_TABLE_Local;    //relevant for local table, each branch has its own fsm.
    vector<uint32_t> FSM_Table_Global;           //relevant for global table.
    uint32_t GlobalHistReg;                      // binary reg saved as its decimal value


    //cont'
    BTB_t(unsigned btbSize, unsigned historySize, unsigned tagSize, unsigned fsmState,
          bool isGlobalHist, bool isGlobalTable, uint32_t Shared) : entries(btbSize), BHR_size(historySize),
                                                                    tag_size(tagSize), start_state(fsmState),
                                                                    isGlobalHist(isGlobalHist),
                                                                    isGlobalTable(isGlobalTable), Shared(Shared),
                                                                    GlobalHistReg(0) {}

    ~BTB_t() = default;
} *BTB;

//Global Variables
BTB our_BTB; //
SIM_stats *our_stats;

//--------------------------------------------Aux functions-----------------------------------------------
uint32_t convertBinaryToDecimal(vector<uint32_t> arr) {
    uint32_t len = arr.size();
    uint32_t res = 0;
    for (uint32_t i = 0; i < len; ++i) {
        if (arr[i]) {
            res += pow(2, i);
        }
    }
    return res;
}

vector<uint32_t> convertDecimalToBinary(uint32_t Dec) {
    vector<uint32_t> res = vector<uint32_t>();
    for (uint32_t i = 0; i < 32; ++i) {
        res.push_back(Dec % 2);
        Dec = Dec / 2;
    }
    return res;
}

bool predict_GHR_gFSM(uint32_t table_index, uint32_t given_tag, uint32_t xor_arg, uint32_t pc,
                      uint32_t *dst) {                           //least space required.

    if (our_BTB->HistoryTable[table_index][0] &&
        (given_tag == our_BTB->HistoryTable[table_index][1])) //if the branch is found in the table and matches the given tag.
    {
        uint32_t fsm_row = 0;
        if (our_BTB->Shared == NOT_USING_SHARE) {
            fsm_row = our_BTB->GlobalHistReg;
        } else {
            fsm_row = (xor_arg) ^ (our_BTB->GlobalHistReg);
        }
        if (our_BTB->FSM_Table_Global[fsm_row] == WT ||
            our_BTB->FSM_Table_Global[fsm_row] == ST) { // predict taken case
            *dst = our_BTB->HistoryTable[table_index][3];
            return true;
        }
    }
    //predict NT case or not found in the table.
    *dst = pc + 4;
    return false;
}

bool predict_LHR_gFSM(uint32_t table_index, uint32_t given_tag, uint32_t xor_arg, uint32_t pc, uint32_t *dst) {
    if (our_BTB->HistoryTable[table_index][0] &&
        (given_tag == our_BTB->HistoryTable[table_index][1])) //if the branch is found in the table and matches the given tag.
    {
        uint32_t fsm_row = 0;
        if (our_BTB->Shared == NOT_USING_SHARE) {
            fsm_row = our_BTB->HistoryTable[table_index][2];
        } else {
            fsm_row = (xor_arg) ^ (our_BTB->HistoryTable[table_index][2]);
        }
        if (our_BTB->FSM_Table_Global[fsm_row] == WT ||
            our_BTB->FSM_Table_Global[fsm_row] == ST) {     // predict taken case
            *dst = our_BTB->HistoryTable[table_index][3];
            return true;
        }
    }
    //predict NT case or not found in the table.
    *dst = pc + 4;
    return false;
}

bool predict_GHR_lFSM(uint32_t table_index, uint32_t given_tag, uint32_t pc, uint32_t *dst) {
    if (our_BTB->HistoryTable[table_index][0] &&
        (given_tag == our_BTB->HistoryTable[table_index][1])) //if the branch is found in the table and matches the given tag.
    {
        if (our_BTB->FSM_TABLE_Local[table_index][our_BTB->GlobalHistReg] == ST ||
            our_BTB->FSM_TABLE_Local[table_index][our_BTB->GlobalHistReg] == WT) {
            *dst = our_BTB->HistoryTable[table_index][3];
            return true;
        }
    }
    //predict NT case or not found in the table.
    *dst = pc + 4;
    return false;
}

bool predict_LHR_lFSM(uint32_t table_index, uint32_t given_tag, uint32_t pc, uint32_t *dst) {
    if (our_BTB->HistoryTable[table_index][0] &&
        (given_tag == our_BTB->HistoryTable[table_index][1])) //if the branch is found in the table and matches the given tag
    {
        if (our_BTB->FSM_TABLE_Local[table_index][our_BTB->HistoryTable[table_index][2]] == ST ||
            our_BTB->FSM_TABLE_Local[table_index][our_BTB->HistoryTable[table_index][2]] == WT) {
            *dst = our_BTB->HistoryTable[table_index][3];
            return true;
        }
    }
    //predict NT case or not found in the table.
    *dst = pc + 4;
    return false;
}


uint32_t updateFsmState(uint32_t state, bool taken) {

    switch (state) {
        case SNT: {
            return taken ? WNT : SNT;
        }
        case WNT: {
            return taken ? WT : SNT;
        }
        case WT: {
            return taken ? ST : WNT;
        }
        case ST: {
            return taken ? ST : WT;
        }
        default:
            return state;
    }
}

uint32_t updateHistory(uint32_t history, bool taken) {
    unsigned int mask = pow(2, our_BTB->BHR_size) - 1;      //mask = 111...11 --> the number of the bits = size of the history register.
    history = (history << 1) & mask;
    if (taken == 1) {
        history |= 1;
    }
    return history;
}


//--------------------------------------------------------------------------------------------------------
int BP_init(unsigned btbSize, unsigned historySize, unsigned tagSize, unsigned fsmState,
            bool isGlobalHist, bool isGlobalTable, int Shared) {

    our_BTB = new BTB_t(btbSize, historySize, tagSize, fsmState, isGlobalHist, isGlobalTable, Shared);

    //history table and fsm tables initialized according to local / global cases.
    our_BTB->HistoryTable = vector<vector<uint32_t>>(btbSize,vector<uint32_t>(MAIN_TABLE_COL,0)); //needed in both local and global history.

    if (isGlobalTable) {
        our_BTB->FSM_Table_Global = vector<uint32_t>(pow(2, historySize), fsmState);
    } else { //local table.
        our_BTB->FSM_TABLE_Local = vector<vector<uint32_t>>(btbSize, vector<uint32_t>(pow(2, historySize), fsmState));
    }


    our_stats = new SIM_stats;
    our_stats->br_num = 0;
    our_stats->flush_num = 0;
    our_stats->size = 0;

    // modify stats struct field sizes accordingly
    if ((our_BTB->isGlobalTable) && (our_BTB->isGlobalHist)) {
        our_stats->size =
                our_BTB->entries * (our_BTB->tag_size + VALID_BIT_S + TARGET_PC) + 2 * pow(2, our_BTB->BHR_size) +
                our_BTB->BHR_size;
    }
    if (!(our_BTB->isGlobalTable) && (our_BTB->isGlobalHist)) {
        our_stats->size =
                our_BTB->entries * (our_BTB->tag_size + VALID_BIT_S + TARGET_PC + 2 * pow(2, our_BTB->BHR_size)) +
                our_BTB->BHR_size;
    }
    if (our_BTB->isGlobalTable && !(our_BTB->isGlobalHist)) {
        our_stats->size = our_BTB->entries * (our_BTB->tag_size + VALID_BIT_S + TARGET_PC + our_BTB->BHR_size) +
                          2 * pow(2, our_BTB->BHR_size);
    }
    if (!(our_BTB->isGlobalTable) && !(our_BTB->isGlobalHist)) {
        our_stats->size = our_BTB->entries * (our_BTB->tag_size + VALID_BIT_S + TARGET_PC + our_BTB->BHR_size +
                                              2 * pow(2, our_BTB->BHR_size));
    }
    return 0;
}



bool BP_predict(uint32_t pc, uint32_t *dst) {

    vector<uint32_t> tmp_array(32, 0);          //represent the pc as a binary number by arranging the bits in a vector.
    tmp_array = convertDecimalToBinary(pc);
    uint32_t x = (uint32_t) log2(our_BTB->entries); //x = number of bits needed to define an entry in the BTB.

    //splitting the vector of the pc
    vector<uint32_t> share_lsb_binary(tmp_array.begin() + 2, tmp_array.begin() + 2 + our_BTB->BHR_size);
    vector<uint32_t> share_mid_binary(tmp_array.begin() + 16, tmp_array.begin() + 16 + our_BTB->BHR_size);
    vector<uint32_t> table_row_binary(tmp_array.begin() + 2, tmp_array.begin() + 2 + x);


    uint32_t end_of_tag = 0;
    if (our_BTB->tag_size <= 30 - x) {                  //the tag does not overstep the pc size (32 bits).
        end_of_tag = 30 - x - our_BTB->tag_size;        //number of the last bit of the tag.
    }
    vector<uint32_t> given_tag_binary(tmp_array.begin() + 2 + x, tmp_array.begin() + 2 + (30 - end_of_tag));

    //converting the binary vectors to decimal numbers.
    uint32_t table_row = convertBinaryToDecimal(table_row_binary);
    uint32_t given_tag = convertBinaryToDecimal(given_tag_binary);
    uint32_t lsb_xor = convertBinaryToDecimal(share_lsb_binary);
    uint32_t mid_xor = convertBinaryToDecimal(share_mid_binary);


    //Global History + Global FSM Table
    if ((our_BTB->isGlobalTable) && (our_BTB->isGlobalHist)) {
        if (our_BTB->Shared == USING_SHARE_LSB) {
            return predict_GHR_gFSM(table_row, given_tag, lsb_xor, pc, dst);
        } else return predict_GHR_gFSM(table_row, given_tag, mid_xor, pc, dst);
    }

    //Global History + Local FSM Table
    if (!(our_BTB->isGlobalTable) && (our_BTB->isGlobalHist)) {
        return predict_GHR_lFSM(table_row, given_tag, pc, dst);
    }
    //Local History + Global FSM Table
    if (our_BTB->isGlobalTable && !(our_BTB->isGlobalHist)) {
        if (our_BTB->Shared == USING_SHARE_LSB) {
            return predict_LHR_gFSM(table_row, given_tag, lsb_xor, pc, dst);
        } else return predict_LHR_gFSM(table_row, given_tag, mid_xor, pc, dst);
    }
    //Local History + Local FSM Table
    if (!(our_BTB->isGlobalTable) && !(our_BTB->isGlobalHist)) {
        return predict_LHR_lFSM(table_row, given_tag, pc, dst);
    }

    return false;
}



void BP_update(uint32_t pc, uint32_t targetPc, bool taken, uint32_t pred_dst) {
    our_stats->br_num++;
    vector<uint32_t> tmp_array(32, 0);              //represent the pc as a binary number by arranging the bits in a vector.
    tmp_array = convertDecimalToBinary(pc);
    uint32_t x = (uint32_t) log2(our_BTB->entries);         //x = number of bits needed to define an entry in the BTB.

    //splitting the vector of the pc
    vector<uint32_t> share_lsb_binary(tmp_array.begin() + 2, tmp_array.begin() + 2 + our_BTB->BHR_size);
    vector<uint32_t> share_mid_binary(tmp_array.begin() + 16, tmp_array.begin() + 16 + our_BTB->BHR_size);
    vector<uint32_t> table_row_binary(tmp_array.begin() + 2, tmp_array.begin() + 2 + x);


    uint32_t end_of_tag = 0;
    if (our_BTB->tag_size <= 30 - x) {                     //the tag does not overstep the pc size (32 bits).
        end_of_tag = 30 - x - our_BTB->tag_size;           //number of the last bit of the tag.
    }
    vector<uint32_t> given_tag_binary(tmp_array.begin() + 2 + x, tmp_array.begin() + 2 + (30 - end_of_tag));

    //converting the binary vectors to decimal numbers.
    uint32_t table_row = convertBinaryToDecimal(table_row_binary);
    uint32_t given_tag = convertBinaryToDecimal(given_tag_binary);
    uint32_t lsb_xor = convertBinaryToDecimal(share_lsb_binary);
    uint32_t mid_xor = convertBinaryToDecimal(share_mid_binary);



    /*----------- Local history & Local fsm--------------------*/
    if (!(our_BTB->isGlobalHist) && !(our_BTB->isGlobalTable)) {
        if (our_BTB->HistoryTable[table_row][0] &&
            (our_BTB->HistoryTable[table_row][1] == given_tag)) {                               // if branch exists in the table

            uint32_t curr_hist = our_BTB->HistoryTable[table_row][2];
            our_BTB->FSM_TABLE_Local[table_row][curr_hist] = updateFsmState(our_BTB->FSM_TABLE_Local[table_row][curr_hist], taken);                     // update the state in matching local fsm
            our_BTB->HistoryTable[table_row][2] = updateHistory(curr_hist, taken);
        } else {                                                                                //branch doesn't exist
            uint32_t new_history = 0;
            for (uint32_t i = 0; i < our_BTB->FSM_TABLE_Local[table_row].size(); ++i) {
                our_BTB->FSM_TABLE_Local[table_row][i] = our_BTB->start_state;                  //initializing the local fsm
            }
            our_BTB->FSM_TABLE_Local[table_row][new_history] = updateFsmState(our_BTB->start_state, taken);
            our_BTB->HistoryTable[table_row][2] = taken;
        }

    }


    /*----------- Global history & Local fsm--------------------*/
    if ((our_BTB->isGlobalHist) && !(our_BTB->isGlobalTable)) {
        if (our_BTB->HistoryTable[table_row][0] &&
            (our_BTB->HistoryTable[table_row][1] == given_tag)) {                                   // if branch exists in the table

            uint32_t curr_history = our_BTB->GlobalHistReg;
            our_BTB->FSM_TABLE_Local[table_row][curr_history] = updateFsmState(our_BTB->FSM_TABLE_Local[table_row][curr_history], taken);
        } else {                                                                                    //branch doesn't exist
            for (uint32_t i = 0; i < our_BTB->FSM_TABLE_Local[table_row].size(); ++i) {
                our_BTB->FSM_TABLE_Local[table_row][i] = our_BTB->start_state;                       //initializing the local fsm
            }
            our_BTB->FSM_TABLE_Local[table_row][our_BTB->GlobalHistReg] = updateFsmState(our_BTB->start_state, taken);
        }
        our_BTB->GlobalHistReg = updateHistory(our_BTB->GlobalHistReg, taken);
    }



    /*----------- Global history & Global fsm--------------------*/
    if ((our_BTB->isGlobalHist) && (our_BTB->isGlobalTable)) {
        uint32_t curr_hist = our_BTB->GlobalHistReg;
        if (our_BTB->Shared == USING_SHARE_LSB) {
            uint32_t tmp = 0;
            tmp = curr_hist;
            curr_hist = (tmp ^ lsb_xor);
        } else if (our_BTB->Shared == USING_SHARE_MID) {
            uint32_t tmp = 0;
            tmp = curr_hist;
            curr_hist = (tmp ^ mid_xor);
        } else {    //NOT_USING_SHARE
            curr_hist = our_BTB->GlobalHistReg;
        }
        our_BTB->FSM_Table_Global[curr_hist] = updateFsmState(our_BTB->FSM_Table_Global[curr_hist], taken);
        our_BTB->GlobalHistReg = updateHistory(our_BTB->GlobalHistReg, taken);
    }


    /*----------- Local history & Global fsm--------------------*/
    if (!(our_BTB->isGlobalHist) && (our_BTB->isGlobalTable)) {
        uint32_t curr_history = 0;
        if (our_BTB->HistoryTable[table_row][0] &&
            (our_BTB->HistoryTable[table_row][1] == given_tag)) {                                // if branch exists in the table
            curr_history = our_BTB->HistoryTable[table_row][2];
            our_BTB->HistoryTable[table_row][2] = updateHistory(curr_history, taken);
        } else {                                                                                 //branch doesn't exist
            our_BTB->HistoryTable[table_row][2] = updateHistory(curr_history, taken);
        }
        if (our_BTB->Shared == USING_SHARE_LSB) {
            uint32_t tmp = 0;
            tmp = curr_history;
            curr_history = (tmp ^ lsb_xor);
        }
        if (our_BTB->Shared == USING_SHARE_MID) {
            uint32_t tmp = 0;
            tmp = curr_history;
            curr_history = (tmp ^ mid_xor);
        }
        our_BTB->FSM_Table_Global[curr_history] = updateFsmState(our_BTB->FSM_Table_Global[curr_history], taken);
    }

    //update the history table
    our_BTB->HistoryTable[table_row][0] = 1;
    our_BTB->HistoryTable[table_row][1] = given_tag;
    our_BTB->HistoryTable[table_row][3] = targetPc;

    if ((taken && (targetPc != pred_dst)) || (!taken && (pred_dst != pc + 4)))
        our_stats->flush_num++;

    return;

}



void BP_GetStats(SIM_stats *curStats) {
    curStats->size = our_stats->size;
    curStats->flush_num = our_stats->flush_num;
    curStats->br_num = our_stats->br_num;
    delete our_BTB;
    delete our_stats;
    return;
}




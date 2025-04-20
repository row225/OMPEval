#include "omp/EquityCalculator.h"
#include "omp/CardRange.h"
#include <iostream>
#include <vector>
#include <string>
#include <sstream>

using namespace std;
using namespace omp;

// Helper to join board cards into a comma-separated string
string joinBoard(const vector<string>& cards) {
    stringstream ss;
    for (size_t i = 0; i < cards.size(); ++i) {
        if (i > 0) ss << ",";
        ss << cards[i];
    }
    return ss.str();
}

// Run equity calculation and print results
void runSimulation(const string& heroHand, const vector<CardRange>& opponentRanges, const vector<string>& board) {
    vector<CardRange> allRanges = {CardRange(heroHand)};
    allRanges.insert(allRanges.end(), opponentRanges.begin(), opponentRanges.end());

    uint64_t boardMask = board.empty() ? 0 : CardRange::getCardMask(joinBoard(board));
    EquityCalculator eq;
    eq.start(allRanges, boardMask, 0, false, 1e-5); // 0.001% margin
    eq.wait();

    auto r = eq.getResults();
    double win = 100.0 * r.equity[0];
    double tie = 100.0 * r.equity[2];
    double loss = 100.0 * (1.0 - r.equity[0] - r.equity[2]);

    cout << "Win: " << win << "%, Tie: " << tie << "%, Loss: " << loss << "%" << endl;
}

int main() {
    string card1, card2;
    cout << "Enter your two hole cards (e.g., As Ks): ";
    cin >> card1 >> card2;
    string heroHand = card1 + card2;

    int numOpponents;
    cout << "Enter number of opponents: ";
    cin >> numOpponents;

    vector<CardRange> opponentRanges;
    cin.ignore(); // Clear newline
    for (int i = 0; i < numOpponents; ++i) {
        string range;
        cout << "Enter range for Opponent " << (i + 1) << " (e.g., QQ+, AKs, random): ";
        getline(cin, range);
        opponentRanges.push_back(CardRange(range));
    }

    vector<string> board;

    cout << "\n[Step 1: Pre-Flop]" << endl;
    runSimulation(heroHand, opponentRanges, board);

    // Flop
    cout << "\n[Step 2: Flop]" << endl;
    cout << "Enter 3 flop cards (e.g., 2h 5h Qc): ";
    for (int i = 0; i < 3; ++i) {
        string card;
        cin >> card;
        board.push_back(card);
    }
    runSimulation(heroHand, opponentRanges, board);

    // Turn
    cout << "\n[Step 3: Turn]" << endl;
    cout << "Enter turn card: ";
    string turn;
    cin >> turn;
    board.push_back(turn);
    runSimulation(heroHand, opponentRanges, board);

    // River
    cout << "\n[Step 4: River]" << endl;
    cout << "Enter river card: ";
    string river;
    cin >> river;
    board.push_back(river);

    cout << "\nFinal Result → ";
    runSimulation(heroHand, opponentRanges, board);

    return 0;
}

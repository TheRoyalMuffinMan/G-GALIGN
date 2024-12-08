#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <mutex>
#include <thread>
#include "include/cmd.hpp"
#include "include/fasta.hpp"
#include "include/globals.hpp"

// Set to 0 to disable debugging
#define DDEBUG 1
#define nTHREADS 3 // CHANGE TO INCREASE MAX NUM THREADS / Matrix Sub-Blocks

// GLOBAL VALUES //
// Global 2D Array For Alignment Scores & Their Directions & ACSII Symbol
int64_t** matrixValues;
char** matrixDir;
char** matrixSym;

// Alignment Scoring Values & Sequences
int gap_penalty;
int mismatch;
int match;
int queryL;
int refL;
std::string rawQuery;
std::string rawReference;
bool ignoreSEGaps = false; //args.ignore_outer_gaps;

// Mutex & Semaphores 
std::mutex exclusive; // Mutex Lock for Alignment
std::mutex threadLock; // Lock for adding threads to the thread list/pool

// Thread List/Pool
std::thread threadList[nTHREADS];

// Thread Block Class
class threadBlock {
    public:
        int startRow = -1;
        int startCol = -1;
        int endCol = -1;
        int endRow = -1;
};

// Thread Alignment Function
void blockAlign(int startRow, int stopRow, int startCol, int stopCol)
{
    // Loop over all matrix values starting at (1,1)
    for(int q = startRow; q < stopRow; q++)
    {
        for(int r = startCol; r < stopCol; r++) 
        {
            // Score Variables for each of the 3 possible directions
            int64_t leftTile = 0;
            int64_t upTile = 0;
            int64_t diaTile = 0;

            // LEFT SCORE
            // If we are on the last column & ignoring start/end gap penalties
            if (r == (refL) && ignoreSEGaps == true)
            {
                leftTile = matrixValues[q][r-1] + 0;
            }
            // Otherwise add normal gap penalty to left tile's score
            else
            {
                leftTile = matrixValues[q][r-1] + gap_penalty;
            }

            // ABOVE SCORE
            // If we are on the last row & ignoring start/end gap penalties
            if (q == (queryL) && ignoreSEGaps == true)
            {
                upTile = matrixValues[q-1][r] + 0;
            }

            // Otherwise add normal gap penalty to above tile's score
            else
            {
                upTile = matrixValues[q-1][r] + gap_penalty;
            }

            // DIAGONAL SCORE
            char maxSym = '_';
            if (rawQuery[q-1] == rawReference[r-1]) // Match Num - 1, as the refrences start at 0 to M but the matrix goes from 0 to M+1 due to the additional init row/col
            {
                diaTile = matrixValues[q-1][r-1] + match; // Add match score to upper diagonal's score
                maxSym = '|'; // Match Symbol
            }
            else
            {
                diaTile = matrixValues[q-1][r-1] + mismatch; // Add mismatch score to upper diagonal's score
                maxSym = 'x'; // Mismatch Symbol
            }

            // Save the Best Score
            int64_t maxScore = diaTile; // Default Diagonal Tile
            char maxDir = 'D'; // Direction we used to calculate the score for this tile
            // If the score from coming from above is bigger switch the max score, and maxDir values
            if (upTile > maxScore)
            {
                maxScore = upTile;
                maxDir = 'U';
                maxSym = ' '; // Store symbol for a gap
            }
            // If the score from the left gap addition is bigger
            if (leftTile > maxScore)
            {
                maxScore = leftTile;
                maxDir = 'L';
                maxSym = ' ';
            }

            // Get the mutex lock
            exclusive.lock();

            // Store the max score, and its corresponding direction and visualization symbol
            matrixValues[q][r] = maxScore;
            matrixDir[q][r] = maxDir;
            matrixSym[q][r] = maxSym;

            // Release the mutex lock
            exclusive.unlock();
        }
    }
}

// Main Function
int main(int argc, char *argv[]) {
    CommandLineArgs args(argc, argv);
    if (args.parse() != 0) {
        std::exit(EXIT_FAILURE);
    }
    #if (DDEBUG != 0)
        args.print();
    #endif

    Fasta query_fasta(args.query);
    Fasta reference_fasta(args.reference);
    if (query_fasta.read() != 0 || reference_fasta.read() != 0) {
        std::exit(EXIT_FAILURE);
    }
    #if (DDEBUG != 0)
        query_fasta.print();
        reference_fasta.print();
    #endif

    // Only consider the first entry in each file for global alignment
    Gene query = query_fasta.genes[0];
    Gene reference = reference_fasta.genes[0];

    //write_results(10, args.output, reference.id, reference.sequence, "", query.id, query.sequence);

    // Get the arguments
    std::string output_file = args.output;
    gap_penalty = args.gap_penalty;
    mismatch = args.mismatch_penalty;
    match = args.match_score;
    ignoreSEGaps = false; //args.ignore_outer_gaps;

    // The raw query and reference strings
    rawQuery = query.sequence;
    rawReference = reference.sequence;

    // Get the length of both sequences
    queryL = rawQuery.length() + 1;
    refL = rawReference.length() + 1;

    // Create 2D Array For alignment scores & their directions & ACSII Symbol (| = match, x = missmatch, " " = gap, S = start)
    matrixValues = new int64_t*[queryL];
    matrixDir = new char*[queryL];
    matrixSym = new char*[queryL];

    // Manually Building a dynamically allocated 2D Array
    for (int i = 0; i < queryL; i++)
    {
        matrixValues[i] = new int64_t[refL];
        matrixDir[i] = new char[refL];
        matrixSym[i] = new char[refL];
    }

    // Manually Initilizing 2D Arrays
    for (int i = 0; i < queryL; i++)
    {
        for (int j=0; j < refL; j++)
        {
            matrixValues[i][j] = 0; // for each empty array fill with 0s, one for each column in that row
            matrixSym[i][j] = 'Z'; // Init with Zs
            if(i != 0)
            {
                matrixDir[i][j] = 'U'; // Fill with U for Up for direction initilization
            }
            else
            {
                matrixDir[i][j] = 'L'; // Fill with L for Left for direction initilization
            }
        }
    }
    // Init First Row and Col of Start Gaps (Start Tile)
    matrixValues[0][0] = 0;
    matrixSym[0][0] = 'S';
    matrixDir[0][0] = '-';

    // Start the timing Calculations
    auto start = std::chrono::high_resolution_clock::now();

    // GAP ROW AND COLUMN INITILIZATION

    // Init the first column and first row of the matrix, with increasing gap penalty values
    // Init the First Column
    for (int i = 0; i < queryL; i++) // Rows
    {
        // If no start gap penalty, 0s
        if(ignoreSEGaps == true)
        {
            matrixValues[i][0] = 0;
            matrixSym[i][0] = ' ';
        }
        else
        {
            if (i != 0)
            {
                matrixValues[i][0] = gap_penalty + matrixValues[i-1][0];
                matrixSym[i][0] = ' ';
            }
        }
    }

    // Init the First Row
    for (int j = 0; j < refL; j++) // Cols
    {
        // If no start gap penalty, 0s
        if(ignoreSEGaps)
        {
            matrixValues[0][j] = 0;
            matrixSym[0][j] = ' ';
        }
        else
        {
            if (j != 0)
            {
                matrixValues[0][j] = gap_penalty + matrixValues[0][j-1];
                matrixSym[0][j] = ' ';
            }
        }
    }

    // ################################################### START OF THREAD APPROACH ############################################################
    // Calculate the Alignment Scores
    
    // Divide cols 1-m and rows 1-n into nThread number of sections each. (SKIP 0th row/col bc its already filled out with gap values)
    int numRowsToAdd = queryL / (nTHREADS);
    int numColsToAdd = refL / (nTHREADS);

    // The number of active threads
    int numThreads = 0;

    // List of the existing thread blocks
    threadBlock existingBlocks[nTHREADS];

    //std::cout << queryL << " Ref: " << refL << "\n";

    // Add new threads as we move along the anti diagonal sub blocks, 0,0  down to the max number of threads we will need to make (Biggest Diagonal)
    while (numThreads < nTHREADS)
    {
        // List of new matrix blocks for the next set of threads, & Index
        threadBlock newBlocks[nTHREADS];
        int newBlockIndex = 0;

        // For Each Block From the Previous Itteration, calculate the next two blocks
        for (int n = 0; n < numThreads; n++)
        {
            // ROW //
            // Calculate a new block moving along the rows (col values stay the same)
            int newStartRow = existingBlocks[n].endRow;
            int newStopRow = existingBlocks[n].endRow + numRowsToAdd;

            // Check if this block is out of bounds (error state)
            if (newStartRow > queryL || newStopRow > queryL)
            {
                std::cout << "ERROR IN SUB-BLOCK ROW CALCULATION \n";
                return 0;
            }

            // Check if we have already added this sub-block
            bool newRowBlock = true;
            for (int i = 0; i < numThreads + 1; i++)
            {
                if (newBlocks[i].startRow == newStartRow)
                {
                    newRowBlock = false;
                    break;
                }
            }

            // If it is a new block
            if (newRowBlock == true)
            {
                // Check if we are near the end of the string, if so bump out the end column
                if ((queryL - newStopRow) < (nTHREADS + 1))
                {
                    newStopRow = queryL;
                }

                // Save the Block
                newBlocks[newBlockIndex].startCol = existingBlocks[n].startCol;
                newBlocks[newBlockIndex].endCol = existingBlocks[n].endCol;
                newBlocks[newBlockIndex].startRow = newStartRow;
                newBlocks[newBlockIndex].endRow = newStopRow;
                newBlockIndex++;

                // Create a new Thread
                //std::cout << "R : Itteration " << numThreads << "  New Block: " << existingBlocks[n].startCol << " Col End: " << existingBlocks[n].endCol << " Row Start: " << newStartRow << " Row End: " <<  newStopRow << "\n";
                threadLock.lock(); // Get Lock
                threadList[newBlockIndex-1] = std::thread (blockAlign, newStartRow, newStopRow, existingBlocks[n].startCol, existingBlocks[n].endCol);
                threadLock.unlock(); // Release Lock
            }

            // COL //
            // Calculate a new block moving along the cols (row values stay the same)
            int newStartCol = existingBlocks[n].endCol;
            int newStopCol = existingBlocks[n].endCol + numColsToAdd;

            // Check if this block is out of bounds (error state)
            if (newStartCol > refL || newStopCol > refL)
            {
                std::cout << "ERROR IN SUB-BLOCK ROW CALCULATION \n";
                return 0;
            }

            // Check if we have already added this sub-block
            bool newColBlock = true;
            for (int j = 0; j < numThreads + 1; j++)
            {
                if (newBlocks[j].startCol == newStartCol)
                {
                    newColBlock = false;
                    break;
                }
            }

            // If it is a new block
            if (newColBlock == true)
            {
                // Check if we are near the end of the string, if so bump out the end column
                if ((refL - newStopCol) < (nTHREADS + 1))
                {
                    newStopCol = refL;
                }

                // Save the Block
                newBlocks[newBlockIndex].startCol = newStartCol;
                newBlocks[newBlockIndex].endCol = newStopCol;
                newBlocks[newBlockIndex].startRow = existingBlocks[n].startRow;
                newBlocks[newBlockIndex].endRow = existingBlocks[n].endRow;
                newBlockIndex++;

                // Create a new Thread
                //std::cout << "C : Itteration " << numThreads << "  New Block: " << newStartCol << " Col End: " << newStopCol << " Row Start: " << existingBlocks[n].startRow << " Row Stop: " <<  existingBlocks[n].endRow << "\n";
                threadLock.lock();
                threadList[newBlockIndex-1] = std::thread (blockAlign, existingBlocks[n].startRow, existingBlocks[n].endRow, newStartCol, newStopCol);
                threadLock.unlock();
            }
        }

        // Start the process by calculating the first block, if we have no prev threads
        if (numThreads == 0)
        {
            // Add the block to the list of complete blocks
            newBlocks[0].startCol = 1;
            newBlocks[0].startRow = 1;
            newBlocks[0].endCol = numColsToAdd;
            newBlocks[0].endRow = numRowsToAdd;

            // Create the thread
            //std::cout << " Num Col " << numColsToAdd << " Num Rows " << numRowsToAdd << "\n";
            std::thread t0(blockAlign, 1, numRowsToAdd, 1, numRowsToAdd);
            t0.join(); // Wait for Thread to Complete (Only 1 Thread in this Itteration)
        }

        // Wait for all the threads to complete before the next itteration
        for (int x=0; x < newBlockIndex; x++)
        {
            threadList[x].join();
        }

        // Increase the thread count by 1 
        numThreads++;

        // Save the current thread blocks
        for (int z = 0; z < nTHREADS; z++)
        {
            existingBlocks[z] = newBlocks[z];
        }
    }

    //Loop over the remaining diagonals, shrinking the number of active threads until we calculate the last sub block
    while (numThreads > -1)
    {
        // List of new matrix blocks for the next set of threads, & Index
        threadBlock newBlocks[nTHREADS];
        int newBlockIndex = 0;

        // For Each Block From the Previous Itteration, calculate the next two blocks
        for (int n = 0; n < numThreads; n++)
        {
            // ROW //
            // Calculate a new block moving along the rows (col values stay the same)
            int newStartRow = existingBlocks[n].endRow;
            int newStopRow = existingBlocks[n].endRow + numRowsToAdd;

            // Check if this block is NOT out of bounds (error state)
            if (newStopRow <= queryL)
            {
                // Check if we have already added this sub-block
                bool newRowBlock = true;
                for (int i = 0; i < numThreads + 1; i++)
                {
                    if (newBlocks[i].startRow == newStartRow)
                    {
                        newRowBlock = false;
                        break;
                    }
                }

                // If it is a new block
                if (newRowBlock == true)
                {
                    // Check if we are near the end of the string, if so bump out the end column
                    if ((queryL - newStopRow) < (nTHREADS + 1))
                    {
                        newStopRow = queryL;
                    }

                    // Save the Block
                    newBlocks[newBlockIndex].startCol = existingBlocks[n].startCol;
                    newBlocks[newBlockIndex].endCol = existingBlocks[n].endCol;
                    newBlocks[newBlockIndex].startRow = newStartRow;
                    newBlocks[newBlockIndex].endRow = newStopRow;
                    newBlockIndex++;

                    // Create a new Thread
                    //std::cout << "RD : Itteration " << numThreads << "  New Block: " << existingBlocks[n].startCol << " Col End: " << existingBlocks[n].endCol << " Row Start: " << newStartRow << " Row End: " <<  newStopRow << "\n";
                    threadLock.lock(); // Get Lock
                    threadList[newBlockIndex-1] = std::thread (blockAlign, newStartRow, newStopRow, existingBlocks[n].startCol, existingBlocks[n].endCol);
                    threadLock.unlock(); // Release Lock
                }
            }

            // COL //
            // Calculate a new block moving along the cols (row values stay the same)
            int newStartCol = existingBlocks[n].endCol;
            int newStopCol = existingBlocks[n].endCol + numColsToAdd;

            // Check if this block is NOT out of bounds (error state)
            if (newStopCol <= refL)
            {
                // Check if we have already added this sub-block
                bool newColBlock = true;
                for (int j = 0; j < numThreads + 1; j++)
                {
                    if (newBlocks[j].startCol == newStartCol)
                    {
                        newColBlock = false;
                        break;
                    }
                }

                // If it is a new block
                if (newColBlock == true)
                {
                    // Check if we are near the end of the string, if so bump out the end column
                    if ((refL - newStopCol) < (nTHREADS + 1))
                    {
                        newStopCol = refL;
                    }

                    // Save the Block
                    newBlocks[newBlockIndex].startCol = newStartCol;
                    newBlocks[newBlockIndex].endCol = newStopCol;
                    newBlocks[newBlockIndex].startRow = existingBlocks[n].startRow;
                    newBlocks[newBlockIndex].endRow = existingBlocks[n].endRow;
                    newBlockIndex++;

                    // Create a new Thread
                    //std::cout << "CD : Itteration " << numThreads << "  New Block: " << newStartCol << " Col End: " << newStopCol << " Row Start: " << existingBlocks[n].startRow << " Row Stop: " <<  existingBlocks[n].endRow << "\n";
                    threadLock.lock();
                    threadList[newBlockIndex-1] = std::thread (blockAlign, existingBlocks[n].startRow, existingBlocks[n].endRow, newStartCol, newStopCol);
                    threadLock.unlock();
                }
            }
        }

        // Wait for all the threads to complete before the next itteration
        for (int x=0; x < newBlockIndex; x++)
        {
            threadList[x].join();
        }

        // Decrease the thread count by 1 
        numThreads--;

        // Save the current thread blocks
        for (int z = 0; z < nTHREADS; z++)
        {
            existingBlocks[z] = newBlocks[z];
        }
      
    }

    // End the Timing Calculations
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Execution time: " 
              << duration.count() << " seconds" << std::endl;  // Convert seconds to microseconds

    
    //####################################### ANALYZE RESULTS ########################################
    // Trace the Optimal Path
    // Init Variables for output
    int64_t finalScore = 0;
    std::string refSequence = ""; // Reference Sequence w/gaps
    std::string alignmentVis = "";
    std::string queSequence = ""; // Query Sequence w/gaps

    int indexQ = 0;
    int indexR = 0;

    // If we are ignoring start gaps, different start location
    if(ignoreSEGaps)
    {
        // Starting Point Initially lower right corner
        indexQ = rawQuery.length();
        indexR = rawReference.length();
        int maxScore = matrixValues[indexQ][indexR];

        // Scan last column and last row for a higher score, if so change start indexes to that tile
        // Scan each row of the last column
        for(int q=0; q < rawQuery.length(); q ++)
        {
            int newScore = matrixValues[q][rawReference.length()];
            if (newScore > maxScore)
            {
                maxScore = newScore;
                indexQ = q;
                indexR = rawReference.length();
            }
        }
        // Scan each col of the last row
        for (int r=0; r < rawReference.length(); r++)
        {
            int newScore = matrixValues[rawQuery.length()][r];
            if (newScore > maxScore)
            {
                maxScore = newScore;
                indexQ = rawQuery.length();
                indexR = r;
            }
        }

        // Save the score at the end tile of the path we are traversing (this is the total score of that path)
        finalScore = matrixValues[indexQ][indexR];

        // If we dont start in the lower right corner of the matrix we must add in the necessary end gaps to the 3, strings
        // If we are in the last row but not the last column
        if (indexQ == rawQuery.length() && indexR != rawReference.length())
        {
            int tempR = indexR + 1;
            //int tempQ = indexQ + 1;
            // Traverse the matrix back to the lower right corner
            while (tempR < (rawReference.length() + 1))
            {
                // Add gap to query, and visualization string
                // Add charecters consumed in the reference string
                alignmentVis = alignmentVis + " ";
                refSequence = refSequence + rawReference[tempR - 1];
                queSequence = queSequence + "_";
                tempR = tempR + 1;
            }
        }
        // If we are in the last column but not the last row
        else if(indexR == rawReference.length() && indexQ != rawQuery.length())
        {
            //int tempR = indexR + 1;
            int tempQ = indexQ + 1;
            // Traverse the matrix back to the lower right corner
            while(tempQ < (rawQuery.length() + 1))
            {
                // Add gap to reference, and visualization string
                // Add charecters consumed in the query string
                alignmentVis = alignmentVis + " ";
                queSequence = queSequence + rawQuery[tempQ - 1];
                refSequence = refSequence + "_";
                tempQ = tempQ + 1;
            }
        }
    }
    else
    {
        // Starting Point w/start and end gap penalties
        // Lower right corner of the matrix
        indexQ = rawQuery.length();
        indexR = rawReference.length();

        // Save the score of the path we are going to traverse in reverse
        finalScore = matrixValues[indexQ][indexR];
    }

    // Loop until we reach the start of the matrix (0,0)
    while (indexQ != 0 || indexR !=0)
    {
        char direction = matrixDir[indexQ][indexR]; // Get the direction we came from to get to the current tile
        char symbolM = matrixSym[indexQ][indexR]; // Get the symbol of the current tile (match, gap, or mismatch)

        // If we got here by a gap from the left
        if (direction == 'L')
        {
            // Add the consumed reference charecter
            refSequence = rawReference[indexR-1] + refSequence; // String Concatenation
            // Add a gap in the query sequence
            queSequence = "_" + queSequence;
            // Add a gap in the visualization string
            alignmentVis = " " + alignmentVis;

            // Update Index (move to the left)
            indexR = indexR - 1;
        }
        // If we got to this tile from a gap from above
        else if (direction == 'U')
        {
            // Add a gap in the reference sequence
            refSequence = "_" + refSequence; //# String Concatenation
            // Add the consumed query sequence
            queSequence = rawQuery[indexQ-1] + queSequence;
            // Add a gap in the visualization string
            alignmentVis = " " + alignmentVis;

            // Update Index (move up)
            indexQ = indexQ - 1;
        }
        // If we got to this tile from the diagonal
        else if(direction == 'D')
        {
            // Add the consumed charecters in both the reference and query to their strings
            refSequence = rawReference[indexR-1] + refSequence; // String Concatenation
            queSequence = rawQuery[indexQ-1] + queSequence;
            // Add the symbol (match, or mismatch) to the alignment visualization string
            alignmentVis = symbolM + alignmentVis;

            // Update Index (move up and left)
            indexQ = indexQ - 1;
            indexR = indexR - 1;        
        }
    }

    //####################################### WRITE OUTPUT ###########################################
    // Open Output File
    std::ofstream outputF(args.output);

    // Write Overall Score & Alignment Visualization to Output File
    outputF << finalScore << "\n";
    outputF << reference.id << "\n";
    outputF << refSequence << "\n";
    outputF << alignmentVis << "\n";
    outputF << queSequence << "\n";
    outputF << query.id << "\n";

    // Close Output File
    outputF.close();

    // FREE ALLOCATED 2D ARRAYS
    for (int i = 0; i < queryL; i++)
    {
        delete[] matrixValues[i];
        delete[] matrixDir[i];
        delete[] matrixSym[i];
    }

    delete[] matrixValues;
    delete[] matrixDir;
    delete[] matrixSym;

}
function updatedMatrix2 = removeIdenticalRow(matrix1, matrix2, rowIndex)
    % Check if the row index is valid
    if rowIndex > size(matrix1, 1)
        error('The specified row index exceeds the number of rows in matrix1.');
    end

    % Get the specified row from matrix1
    rowToCheck = matrix1(rowIndex, :);
    
    % Initialize the updatedMatrix2 as matrix2
    updatedMatrix2 = matrix2;
    
    % Check if the identical row exists in matrix2
    identicalRowIndex = find(ismember(matrix2, rowToCheck, 'rows'));
    
    % If identical row is found, remove it from matrix2
    if ~isempty(identicalRowIndex)
        updatedMatrix2(identicalRowIndex, :) = [];
    end
end
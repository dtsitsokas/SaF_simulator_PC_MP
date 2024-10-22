function B = calculate_means(A, x)
    % Function to calculate the mean of x consecutive elements of every row of matrix A
    % Input: A - input matrix (mxn)
    %        x - number of elements to take the mean of
    % Output: B - output matrix where each element of every row is the mean 
    % of x consecutive elements of the equivalent row of A


    % Determine the number of segments
    n_segments = floor(length(A(1,:)) / x);

    % Preallocate B for performance
    B = zeros(size(A,1), n_segments);

    % Calculate the means
    for j=1:size(A,1)
        for i = 1:n_segments
            start_idx = (i-1)*x + 1;
            end_idx = i*x;
            B(j,i) = mean(A(j,start_idx:end_idx));
        end
    end
end


function Cstar = tumorOutcome(t, C)

    % Smooth slightly to avoid numerical wiggles
    C_smooth = smoothdata(C, 'movmean', 5);

    % Find local minima
    [mins, ~] = findpeaks(-C_smooth);

    if ~isempty(mins)
        % True local minima exist
        Cstar = min(C_smooth);
    else
        % No local minima → use endpoint
        Cstar = C_smooth(end);
    end
end

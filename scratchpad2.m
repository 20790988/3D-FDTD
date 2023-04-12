clear

out(1) = parfeval(@func,1,2);

wait(out(1));
out.OutputArguments{1}

function out = func(seconds)
    pause(seconds);
    out = seconds;
end
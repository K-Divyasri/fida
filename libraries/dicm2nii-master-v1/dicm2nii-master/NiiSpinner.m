classdef NiiSpinner < matlab.mixin.SetGet
% matlab.mixin.SetGet (not plain handle) so generic get(obj,'Prop') /
% set(obj,'Prop',val) function-call syntax works and routes through the
% property accessors below -- nii_viewer.m calls spinners both ways
% (obj.Value and get(obj,'Value')/set(obj,'Value',v)).
% Native (pure-uicontrol) replacement for the MJSpinner/javacomponent-based
% number spinner nii_viewer.m used to build via the old java_spinner() local
% function. Exposes the same property/method surface (Value, Model,
% setValue/getValue, setEnabled, setFont, setEditor, StateChangedCallback,
% ToolTipText) so existing call sites elsewhere in nii_viewer.m are unchanged.
%
% h = NiiSpinner(pos, val, parent, fmt, helpTxt)
%   pos: [left bottom width height]
%   val: [curVal min max step]
%   fmt: '#' for integer, or '#.#', '#.##', ... for decimals

    properties
        Minimum = -Inf
        Maximum = Inf
        StepSize = 1
        StateChangedCallback = ''
    end
    properties (Dependent)
        Value
        Model
        Enable
        ToolTipText
    end
    properties (Access = private)
        Value_ = 0
        NumFmt = '#'
        hEdit
        hUp
        hDown
    end

    methods
        function obj = NiiSpinner(pos, val, parent, fmt, helpTxt)
            if nargin < 4 || isempty(fmt), fmt = '#'; end
            if nargin < 5, helpTxt = ''; end
            obj.NumFmt = fmt;
            obj.Minimum = val(2);
            obj.Maximum = val(3);
            obj.StepSize = val(4);

            x = pos(1); y = pos(2); w = pos(3); h = pos(4);
            bw = 14; % width of the tiny up/down buttons
            obj.hEdit = uicontrol(parent, 'Style', 'edit', ...
                'Position', [x y w-bw h], 'BackgroundColor', 'w', ...
                'HorizontalAlignment', 'center', 'FontSize', 8, ...
                'TooltipString', helpTxt, 'Callback', @(s,e)obj.editCallback());
            obj.hUp = uicontrol(parent, 'Style', 'pushbutton', 'FontSize', 6, ...
                'String', char(9650), 'Position', [x+w-bw y+ceil(h/2) bw floor(h/2)], ...
                'TooltipString', helpTxt, 'Callback', @(s,e)obj.stepBy(1));
            obj.hDown = uicontrol(parent, 'Style', 'pushbutton', 'FontSize', 6, ...
                'String', char(9660), 'Position', [x+w-bw y bw ceil(h/2)], ...
                'TooltipString', helpTxt, 'Callback', @(s,e)obj.stepBy(-1));

            obj.Value_ = val(1);
            obj.updateDisplay();
        end

        function v = get.Value(obj), v = obj.Value_; end
        function set.Value(obj, v)
            v = max(obj.Minimum, min(obj.Maximum, v));
            obj.Value_ = v;
            obj.updateDisplay();
            obj.fireCallback();
        end

        function m = get.Model(obj), m = obj; end % Model.Maximum/StepSize alias self

        function s = get.Enable(obj)
            s = get(obj.hEdit, 'Enable');
        end
        function set.Enable(obj, tf)
            if islogical(tf) || isnumeric(tf)
                s = 'off'; if tf, s = 'on'; end
            else
                s = char(tf);
            end
            set([obj.hEdit obj.hUp obj.hDown], 'Enable', s);
        end

        function s = get.ToolTipText(obj)
            s = get(obj.hEdit, 'TooltipString');
        end
        function set.ToolTipText(obj, s)
            set([obj.hEdit obj.hUp obj.hDown], 'TooltipString', s);
        end

        function setValue(obj, v), obj.Value = v; end
        function v = getValue(obj), v = obj.Value_; end
        function setEnabled(obj, tf), obj.Enable = tf; end
        function setFont(obj, varargin) %#ok<INUSD> % native default font is fine
        end
        function setEditor(obj, fmt) % was: javaObject('javax.swing.JSpinner$NumberEditor',...)
            obj.NumFmt = fmt;
            obj.updateDisplay();
        end
    end

    methods (Access = private)
        function stepBy(obj, dir)
            obj.Value = obj.Value_ + dir * obj.StepSize;
        end
        function editCallback(obj)
            v = str2double(get(obj.hEdit, 'String'));
            if isnan(v), obj.updateDisplay(); return; end % revert invalid entry
            obj.Value = v;
        end
        function updateDisplay(obj)
            if isempty(obj.hEdit) || ~isvalid(obj.hEdit), return; end
            set(obj.hEdit, 'String', obj.formatValue(obj.Value_));
        end
        function s = formatValue(obj, v)
            dotPos = strfind(obj.NumFmt, '.');
            if isempty(dotPos)
                s = sprintf('%d', round(v));
            else
                nDec = numel(obj.NumFmt) - dotPos;
                s = sprintf(['%.' num2str(nDec) 'f'], v);
            end
        end
        function fireCallback(obj)
            cb = obj.StateChangedCallback;
            if isempty(cb), return; end
            if iscell(cb), feval(cb{1}, obj, [], cb{2:end});
            else, feval(cb, obj, []);
            end
        end
    end
end

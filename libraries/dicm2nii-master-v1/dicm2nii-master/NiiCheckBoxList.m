classdef NiiCheckBoxList < handle
% Native (pure-uicontrol) replacement for the JIDE CheckBoxList +
% MJScrollPane / javacomponent combo nii_viewer.m used for its file/overlay
% list. Each row is an independent checkbox (shown/hidden toggle) next to a
% clickable label (selects the row for editing its display parameters),
% matching the original's "click checkbox to toggle, click row to select"
% behavior. Exposes the same method/property names the JIDE object was
% called with (getSelectedIndex, getModel, getCheckBoxListSelectedIndices,
% etc.) so the ~40 call sites elsewhere in nii_viewer.m are unchanged.
%
% obj = NiiCheckBoxList(parent, pos, itemName, tooltipTxt, focusTarget)

    properties
        ValueChangedCallback = ''  % fires on setSelectedIndex / row click ("file" cmd)
        ToggleCallback = ''        % fires on single-index check/uncheck ("toggle" cmd)
        ModelChangedCallback = ''  % fires on any add/remove/rename ("width" cmd)
    end
    properties (SetAccess = private)
        SelectionBackground = [0.31 0.51 0.74] % plain RGB triple (was a java.awt.Color)
        Panel % the uipanel backing this widget; callers resize it directly (was hs.scroll)
    end
    properties (Access = private)
        Items = {}
        Checked = false(1,0)
        SelIndex = -1   % 0-based, -1 = none
        AnchorIdx = -1  % 0-based index of the last checkbox toggled
        RowH = 16
        ScrollOffset = 0
        RowHandles = []
        Slider
        Tooltip = ''
        FocusTarget
    end

    methods
        function obj = NiiCheckBoxList(parent, pos, itemName, tooltipTxt, focusTarget)
            obj.Panel = uipanel(parent, 'Units', 'pixels', 'Position', pos, ...
                'BorderType', 'none', 'Tag', 'NiiCheckBoxList');
            if nargin > 3, obj.Tooltip = tooltipTxt; end
            if nargin > 4, obj.FocusTarget = focusTarget; end
            obj.Items = {itemName};
            obj.Checked = true;
            obj.SelIndex = 0;
            obj.redraw();
        end

        % ---- selection ----
        function idx = getSelectedIndex(obj)
            idx = obj.SelIndex;
        end
        function setSelectedIndex(obj, idx0)
            n = numel(obj.Items);
            if idx0 < 0 || idx0 > n-1, idx0 = -1; end
            if idx0 == obj.SelIndex, return; end
            obj.SelIndex = idx0;
            obj.redraw();
            obj.invoke(obj.ValueChangedCallback);
        end

        % ---- model (list-model surface; getModel() returns self since
        % this class already owns list-model state directly) ----
        function m = getModel(obj)
            m = obj;
        end
        function n = size(obj)
            n = numel(obj.Items);
        end
        function s = get(obj, idx0)
            s = obj.Items{idx0+1};
        end
        function set(obj, idx0, str)
            obj.Items{idx0+1} = str;
            obj.redraw();
            obj.invoke(obj.ModelChangedCallback);
        end
        function remove(obj, idx0)
            n = numel(obj.Items);
            if idx0 < 0 || idx0 >= n, return; end
            obj.Items(idx0+1) = [];
            obj.Checked(idx0+1) = [];
            if obj.SelIndex == idx0, obj.SelIndex = -1;
            elseif obj.SelIndex > idx0, obj.SelIndex = obj.SelIndex - 1;
            end
            obj.redraw();
            obj.invoke(obj.ModelChangedCallback);
        end
        function insertElementAt(obj, str, idx0)
            n = numel(obj.Items);
            idx0 = max(0, min(idx0, n));
            obj.Items = [obj.Items(1:idx0) {str} obj.Items(idx0+1:end)];
            obj.Checked = [obj.Checked(1:idx0) false obj.Checked(idx0+1:end)];
            if obj.SelIndex >= idx0, obj.SelIndex = obj.SelIndex + 1; end
            obj.redraw();
            obj.invoke(obj.ModelChangedCallback);
        end
        function c = toArray(obj)
            c = obj.Items;
        end

        % ---- checkbox (shown/hidden) state ----
        function idx0 = getCheckBoxListSelectedIndices(obj)
            idx0 = find(obj.Checked) - 1;
        end
        function setCheckBoxListSelectedIndices(obj, idx0arr)
            obj.Checked = false(1, numel(obj.Items));
            idx0arr = idx0arr(idx0arr >= 0 & idx0arr < numel(obj.Items));
            obj.Checked(idx0arr+1) = true;
            obj.redraw(); % no callback fired: bulk/programmatic use only
        end
        function addCheckBoxListSelectedIndex(obj, idx0)
            if idx0 < 0 || idx0 >= numel(obj.Items), return; end
            obj.Checked(idx0+1) = true;
            obj.AnchorIdx = idx0;
            obj.redraw();
            obj.invoke(obj.ToggleCallback);
        end
        function removeCheckBoxListSelectedIndex(obj, idx0)
            if idx0 < 0 || idx0 >= numel(obj.Items), return; end
            obj.Checked(idx0+1) = false;
            obj.AnchorIdx = idx0;
            obj.redraw();
            obj.invoke(obj.ToggleCallback);
        end
        function idx0 = getAnchorSelectionIndex(obj)
            idx0 = obj.AnchorIdx;
        end

        % ---- misc ----
        function updateUI(obj)
            obj.redraw();
        end
        function w = getPreferredWidth(obj)
            if isempty(obj.Items), w = 60; return; end
            tmp = uicontrol(obj.Panel, 'Style', 'text', 'Units', 'pixels', ...
                'Visible', 'off', 'FontSize', 8);
            maxW = 0;
            for i = 1:numel(obj.Items)
                set(tmp, 'String', obj.Items{i});
                ext = get(tmp, 'Extent');
                maxW = max(maxW, ext(3));
            end
            delete(tmp);
            w = maxW + 22; % checkbox width + padding
        end
    end

    methods (Access = private)
        function invoke(obj, cb)
            % source arg = obj: the 'toggle' case reads h.getAnchorSelectionIndex
            if isempty(cb), return; end
            if iscell(cb), feval(cb{1}, obj, [], cb{2:end});
            else, feval(cb, obj, []);
            end
        end
        function focusAway(obj)
            if ~isempty(obj.FocusTarget) && isvalid(obj.FocusTarget)
                try, uicontrol(obj.FocusTarget); catch, end %#ok<TRYNC>
            end
        end
        function onCheck(obj, idx0, h)
            obj.Checked(idx0+1) = logical(get(h, 'Value'));
            obj.AnchorIdx = idx0;
            obj.focusAway();
            obj.invoke(obj.ToggleCallback);
        end
        function onSelect(obj, idx0)
            obj.setSelectedIndex(idx0);
            obj.focusAway();
        end
        function onSlide(obj, h)
            n = numel(obj.Items);
            panelPos = get(obj.Panel, 'Position');
            visibleRows = max(1, floor(panelPos(4) / obj.RowH));
            maxOffset = max(0, n - visibleRows);
            obj.ScrollOffset = maxOffset - round(get(h, 'Value'));
            obj.redraw();
        end

        function redraw(obj)
            if isempty(obj.Panel) || ~isvalid(obj.Panel), return; end
            old = obj.RowHandles;
            old = old(isgraphics(old));
            delete(old);
            obj.RowHandles = gobjects(0);

            n = numel(obj.Items);
            panelPos = get(obj.Panel, 'Position');
            w = panelPos(3); ht = panelPos(4);
            visibleRows = max(1, floor(ht / obj.RowH));
            needScroll = n > visibleRows;
            sw = 0; if needScroll, sw = 12; end
            maxOffset = max(0, n - visibleRows);
            obj.ScrollOffset = max(0, min(obj.ScrollOffset, maxOffset));

            if needScroll
                if isempty(obj.Slider) || ~isgraphics(obj.Slider)
                    obj.Slider = uicontrol(obj.Panel, 'Style', 'slider', ...
                        'Callback', @(h,~)obj.onSlide(h));
                end
                set(obj.Slider, 'Position', [w-sw 0 sw ht], ...
                    'Min', 0, 'Max', maxOffset, ...
                    'Value', maxOffset - obj.ScrollOffset, ...
                    'SliderStep', min(1, [1 max(1,visibleRows)] / max(1,maxOffset)), ...
                    'Visible', 'on');
            elseif ~isempty(obj.Slider) && isgraphics(obj.Slider)
                set(obj.Slider, 'Visible', 'off');
            end

            rowIdx0 = obj.ScrollOffset + (0:visibleRows-1);
            rowIdx0 = rowIdx0(rowIdx0 < n);
            labelW = max(1, w - sw - 19);
            for k = 1:numel(rowIdx0)
                i0 = rowIdx0(k);
                y = ht - k*obj.RowH;
                isSel = (i0 == obj.SelIndex);
                bg = get(obj.Panel, 'BackgroundColor');
                if isSel, bg = obj.SelectionBackground; end
                hChk = uicontrol(obj.Panel, 'Style', 'checkbox', ...
                    'Value', obj.Checked(i0+1), 'BackgroundColor', bg, ...
                    'Position', [1 y 16 obj.RowH], ...
                    'Callback', @(h,~)obj.onCheck(i0,h));
                hLbl = uicontrol(obj.Panel, 'Style', 'pushbutton', ...
                    'String', obj.Items{i0+1}, 'FontSize', 8, ...
                    'Position', [18 y labelW obj.RowH], 'BackgroundColor', bg, ...
                    'TooltipString', obj.Tooltip, ...
                    'Callback', @(~,~)obj.onSelect(i0));
                obj.RowHandles(end+1) = hChk; %#ok<AGROW>
                obj.RowHandles(end+1) = hLbl; %#ok<AGROW>
            end
        end
    end
end
